context("One-step and TML estimators produce similar results")
library(data.table)
set.seed(172943)

if (require("sl3")) {
  # Example based on the data-generating mechanism presented in the simulation
  n <- 100
  W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
  A <- rpois(n, lambda = exp(3 + .3 * log(W$W1) - 0.2 * exp(W$W1) * W$W2))
  Y <- rbinom(
    n, 1,
    plogis(-1 + 0.05 * A - 0.02 * A * W$W2 + 0.2 * A * tan(W$W1^2) -
      0.02 * W$W1 * W$W2 + 0.1 * A * W$W1 * W$W2)
  )
  C <- rbinom(n, 1, plogis(rowSums(W) + Y))
  delta_shift <- 2
  EY <- mean(Y)

  # true functional forms
  fitA.0 <- glm(
    A ~ I(log(W1)) + I(exp(W1)):W2,
    family = poisson,
    data = data.frame(A, W)
  )

  fitY.0 <- glm(
    Y ~ A + A:W2 + A:I(tan(W1^2)) + W1:W2 + A:W1:W2,
    family = binomial, data = data.frame(A, W)
  )

  gn.0 <- function(A = A, W = W) {
    dpois(A, lambda = predict(fitA.0, newdata = W, type = "response"))
  }

  Qn.0 <- function(A = A, W = W) {
    predict(
      fitY.0,
      newdata = data.frame(A, W, row.names = NULL),
      type = "response"
    )
  }

  # SL learners to be used for most fits (e.g., IPCW, outcome regression)
  mean_learner <- Lrnr_mean$new()
  glm_learner <- Lrnr_glm$new()
  rf_learner <- Lrnr_ranger$new()
  Q_lib <- Stack$new(mean_learner, glm_learner, rf_learner)
  sl <- Lrnr_sl$new(learners = Q_lib, metalearner = Lrnr_nnls$new())

  # SL learners for fitting the generalized propensity score fit
  hse_learner <- make_learner(Lrnr_density_semiparametric,
    mean_learner = glm_learner
  )
  mvd_learner <- make_learner(Lrnr_density_semiparametric,
    mean_learner = rf_learner,
    var_learner = glm_learner
  )
  g_lib <- Stack$new(hse_learner, mvd_learner)
  cv_selector_density <- Lrnr_cv_selector$new(
    eval_function = loss_loglik_true_cat
  )
  sl_density <- Lrnr_sl$new(
    learners = g_lib,
    metalearner = cv_selector_density
    #metalearner = Lrnr_solnp_density$new()
  )

  # NOTE: using true density like Ivan does
  gn_ext_fitted <- as.data.table(
    lapply(
      c(-delta_shift, 0, delta_shift, 2 * delta_shift),
      function(shift_value) {
        gn_out <- gn.0(A = A + shift_value, W = W)
      }
    )
  )
  setnames(gn_ext_fitted, c("downshift", "noshift", "upshift", "upupshift"))

  # NOTE: should also use true Q for good measure (truth includes interactions)
  Qn_ext_fitted <- as.data.table(
    lapply(c(0, delta_shift), function(shift_value) {
      Qn_out <- Qn.0(A = A + shift_value, W = W)
    })
  )
  setnames(Qn_ext_fitted, c("noshift", "upshift"))

  # fit TMLE
  tmle <- txshift(
    Y = Y, A = A, W = W, delta = delta_shift,
    g_exp_fit_args = list(fit_type = "external"),
    gn_exp_fit_ext = gn_ext_fitted,
    Q_fit = list(fit_type = "external"),
    Qn_fit_ext = Qn_ext_fitted,
    estimator = "tmle"
  )
  tmle_psi <- as.numeric(tmle$psi)

  # fit one-step
  os <- txshift(
    Y = Y, A = A, W = W, delta = delta_shift,
    g_exp_fit_args = list(fit_type = "external"),
    gn_exp_fit_ext = gn_ext_fitted,
    Q_fit_args = list(fit_type = "external"),
    Qn_fit_ext = Qn_ext_fitted,
    estimator = "onestep"
  )
  os_psi <- as.numeric(os$psi)

  # test for reasonable equality between estimators
  test_that("TMLE and one-step implementations match closely", {
    expect_equal(tmle_psi, os_psi, tol = 1e-3)
  })

  # fit TMLE for delta = 0
  tmle_noshift <- txshift(
    Y = Y, A = A, W = W, delta = 0, estimator = "tmle",
    g_exp_fit_args = list(fit_type = "sl", sl_learners_density = sl_density),
    Q_fit_args = list(fit_type = "sl", sl_learners = sl)
  )
  tmle_psi_noshift <- as.numeric(tmle_noshift$psi)

  # fit one-step for delta = 0
  os_noshift <- txshift(
    Y = Y, A = A, W = W, delta = 0, estimator = "onestep",
    g_exp_fit_args = list(fit_type = "sl", sl_learners_density = sl_density),
    Q_fit_args = list(fit_type = "sl", sl_learners = sl)
  )
  os_psi_noshift <- as.numeric(os_noshift$psi)

  # test for reasonable equality between estimators
  test_that("TMLE and one-step match EY very closely for delta = 0", {
    expect_equal(tmle_psi_noshift, EY, tol = 1e-3)
    expect_equal(os_psi_noshift, EY, tol = 1e-3)
  })

  # IPCW-based estimators by adding censoring node
  ipcw_tmle <- txshift(
    W = W, A = A, Y = Y, delta = delta_shift,
    C_samp = C, V = c("W", "Y"),
    estimator = "tmle",
    max_iter = 5,
    samp_fit_args = list(fit_type = "glm"),
    g_exp_fit_args = list(fit_type = "external"),
    gn_exp_fit_ext = gn_ext_fitted[C == 1, ],
    Q_fit_args = list(fit_type = "external"),
    Qn_fit_ext = Qn_ext_fitted[C == 1, ],
    eif_reg_type = "glm"
  )
  ipcw_tmle_psi <- as.numeric(ipcw_tmle$psi)

  ipcw_os <- txshift(
    W = W, A = A, Y = Y, delta = delta_shift,
    C_samp = C, V = c("W", "Y"),
    estimator = "onestep",
    samp_fit_args = list(fit_type = "glm"),
    g_exp_fit_args = list(fit_type = "external"),
    gn_exp_fit_ext = gn_ext_fitted[C == 1, ],
    Q_fit_args = list(fit_type = "external"),
    Qn_fit_ext = Qn_ext_fitted[C == 1, ],
    eif_reg_type = "glm"
  )
  ipcw_os_psi <- as.numeric(ipcw_os$psi)

  # test for reasonable equality between estimators
  test_that("IPCW-augmented TMLE and one-step match reasonably closely", {
    expect_equal(ipcw_tmle_psi, ipcw_os_psi, tol = 1e-3)
  })
}
