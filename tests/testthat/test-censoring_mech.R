context("Estimation with loss to follow-up censoring")
library(data.table)
set.seed(172943)

# Example based on the data-generating mechanism presented in the simulation
n <- 100
W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
A <- rpois(n, lambda = exp(3 + .3 * log(W$W1) - 0.2 * exp(W$W1) * W$W2))
C_cens <- rbinom(n, 1, plogis(rowSums(W) - 1))
Y <- rbinom(
  n, 1,
  plogis(-1 + 0.05 * A - 0.02 * A * W$W2 + 0.2 * A * tan(W$W1^2) -
    0.02 * W$W1 * W$W2 + 0.1 * A * W$W1 * W$W2)
)
C_samp <- rbinom(n, 1, plogis(rowSums(W) + Y))
V <- as.data.table(list(W, Y = Y))
delta_shift <- 2

# fit exposure mechanism with HAL, SL, or GLM
gn_exp_hal <- est_g_exp(
  A = A, W = W,
  delta = delta_shift,
  fit_type = "hal",
  haldensify_args = list(
    n_bins = 5,
    grid_type = "equal_mass",
    lambda_seq = exp(seq(-1, -10, length = 300)),
    # the following arguments are passed to hal9001::fit_hal()
    max_degree = 3, smoothness_orders = 0, num_knots = NULL,
    reduce_basis = 1 / sqrt(length(A))
  )
)

# fit censoring mechanism via SL or GLM
gn_cens_glm <- est_g_cens(
  C_cens = C_cens, A = A, W = W,
  fit_type = "glm"
)

# fit outcome mechanism with GLM
Qn_est_glm <- est_Q(
  Y = Y, C_cens = C_cens, A = A, W = W, delta = delta_shift,
  fit_type = "glm", glm_formula = "Y ~ ."
)

# fit two-phase censoring mechanism with GLM
ipcw_est_glm <- est_samp(V = V, C_samp = C_samp, fit_type = "glm")

# fit TMLE and one-step without SL
tmle_glm <- txshift(
  Y = Y, C_cens = C_cens, A = A, W = W, delta = delta_shift,
  g_exp_fit_args = list(fit_type = "external"),
  gn_exp_fit_ext = gn_exp_hal,
  g_cens_fit_args = list(fit_type = "external"),
  gn_cens_fit_ext = gn_cens_glm,
  Q_fit_args = list(fit_type = "external"),
  Qn_fit_ext = Qn_est_glm,
  estimator = "tmle"
)
os_glm <- txshift(
  Y = Y, C_cens, A = A, W = W, delta = delta_shift,
  g_exp_fit_args = list(fit_type = "external"),
  gn_exp_fit_ext = gn_exp_hal,
  g_cens_fit_args = list(fit_type = "external"),
  gn_cens_fit_ext = gn_cens_glm,
  Q_fit_args = list(fit_type = "external"),
  Qn_fit_ext = Qn_est_glm,
  estimator = "onestep"
)

test_that("TMLE and one-step match for GLM nuisance configurations", {
  tmle_glm_psi <- as.numeric(tmle_glm$psi)
  os_glm_psi <- as.numeric(os_glm$psi)
  expect_equal(tmle_glm_psi, os_glm_psi, tol = 1e-3)
})

# fit IPCW-TMLE and IPCW-one-step with GLM for sampling mechanism estimation
ipcw_tmle_glm <- txshift(
  W = W, A = A, Y = Y, delta = delta_shift,
  C_samp = C_samp, V = c("W", "Y"),
  estimator = "tmle",
  max_iter = 5,
  samp_fit_args = list(fit_type = "external"),
  samp_fit_ext = ipcw_est_glm,
  g_exp_fit_args = list(fit_type = "external"),
  gn_exp_fit_ext = gn_exp_hal[C_samp == 1, ],
  Q_fit_args = list(fit_type = "external"),
  Qn_fit_ext = Qn_est_glm[C_samp == 1, ],
  eif_reg_type = "hal"
)
ipcw_os_glm <- txshift(
  W = W, A = A, Y = Y, delta = delta_shift,
  C_samp = C_samp, V = c("W", "Y"),
  estimator = "onestep",
  samp_fit_args = list(fit_type = "external"),
  samp_fit_ext = ipcw_est_glm,
  g_exp_fit_args = list(fit_type = "external"),
  gn_exp_fit_ext = gn_exp_hal[C_samp == 1, ],
  Q_fit_args = list(fit_type = "external"),
  Qn_fit_ext = Qn_est_glm[C_samp == 1, ],
  eif_reg_type = "hal"
)

test_that("IPCW-TMLE and IPCW-one-step match with GLM for censoring", {
  ipcw_tmle_glm_psi <- as.numeric(ipcw_tmle_glm$psi)
  ipcw_os_glm_psi <- as.numeric(ipcw_os_glm$psi)
  expect_equal(ipcw_tmle_glm_psi, ipcw_os_glm_psi, tol = 5e-2)
})

if (require("sl3")) {
  # SL learners to be used for most fits (e.g., IPCW, outcome regression)
  mean_learner <- Lrnr_mean$new()
  glm_learner <- Lrnr_glm$new()
  rf_learner <- Lrnr_ranger$new()
  Q_lib <- Stack$new(mean_learner, glm_learner, rf_learner)
  sl <- Lrnr_sl$new(learners = Q_lib, metalearner = Lrnr_nnls$new())

  # SL learners for fitting the generalized propensity score fit
  hose_learner <- make_learner(Lrnr_density_semiparametric,
    mean_learner = glm_learner
  )
  hese_learner <- make_learner(Lrnr_density_semiparametric,
    mean_learner = rf_learner,
    var_learner = glm_learner
  )
  g_lib <- Stack$new(hose_learner, hese_learner)
  sl_density <- Lrnr_sl$new(
    learners = g_lib,
    metalearner = Lrnr_solnp_density$new()
  )

  # fit exposure mechanism with SL
  gn_exp_sl <- est_g_exp(
    A = A, W = W,
    delta = delta_shift,
    fit_type = "sl",
    sl_learners_density = sl_density
  )

  # fit censoring mechanism via SL
  gn_cens_sl <- est_g_cens(
    C_cens = C_cens, A = A, W = W,
    fit_type = "sl",
    sl_learners = sl
  )

  # fit outcome mechanism with SL
  Qn_est_sl <- est_Q(
    Y = Y, A = A, W = W, delta = delta_shift,
    fit_type = "sl", sl_learners = sl
  )

  # fit two-phase censoring mechanism with SL
  ipcw_est_sl <- est_samp(
    V = V, C_samp = C_samp, fit_type = "sl",
    sl_learners = sl
  )

  # test for errors when arguments are set inconsistently
  test_that("SL-based nuisance estimation fails without SL library", {
    expect_error(est_g_exp(A = A, W = W, delta = delta_shift, fit_type = "sl"))
    expect_error(est_g_cens(C_cens = C_cens, A = A, W = W, fit_type = "sl"))
    expect_error(est_Q(A = A, W = W, delta = delta_shift, fit_type = "sl"))
    expect_error(est_samp(A = A, W = W, delta = delta_shift, fit_type = "sl"))
  })

  # fit TMLE and one-step for HAL and SL
  tmle_ml <- txshift(
    Y = Y, C_cens = C_cens, A = A, W = W, delta = delta_shift,
    g_exp_fit_args = list(fit_type = "external"),
    gn_exp_fit_ext = gn_exp_sl,
    g_cens_fit_args = list(fit_type = "external"),
    gn_cens_fit_ext = gn_cens_sl,
    Q_fit_args = list(fit_type = "external"),
    Qn_fit_ext = Qn_est_sl,
    estimator = "tmle"
  )
  os_ml <- txshift(
    Y = Y, C_cens = C_cens, A = A, W = W, delta = delta_shift,
    g_exp_fit_args = list(fit_type = "external"),
    gn_exp_fit_ext = gn_exp_sl,
    g_cens_fit_args = list(fit_type = "external"),
    gn_cens_fit_ext = gn_cens_sl,
    Q_fit_args = list(fit_type = "external"),
    Qn_fit_ext = Qn_est_sl,
    estimator = "onestep"
  )

  test_that("TMLE and one-step match for SL nuisance configurations", {
    tmle_ml_psi <- as.numeric(tmle_ml$psi)
    os_ml_psi <- as.numeric(os_ml$psi)
    expect_equal(tmle_ml_psi, os_ml_psi, tol = 1e-1)
  })

  # fit IPCW-TMLE and IPCW-one-step with SL for censoring estimation
  ipcw_tmle_sl <- txshift(
    W = W, A = A, C_cens = C_cens, Y = Y, delta = delta_shift,
    C_samp = C_samp, V = c("W", "Y"),
    estimator = "tmle", max_iter = 5,
    samp_fit_args = list(fit_type = "external"),
    samp_fit_ext = ipcw_est_sl,
    g_exp_fit_args = list(fit_type = "external"),
    gn_exp_fit_ext = gn_exp_sl[C_samp == 1, ],
    g_cens_fit_args = list(fit_type = "external"),
    gn_cens_fit_ext = gn_cens_sl[C_samp == 1],
    Q_fit_args = list(fit_type = "external"),
    Qn_fit_ext = Qn_est_sl[C_samp == 1, ],
    eif_reg_type = "hal"
  )
  ipcw_os_sl <- txshift(
    W = W, A = A, C_cens = C_cens, Y = Y, delta = delta_shift,
    C_samp = C_samp, V = c("W", "Y"),
    estimator = "onestep",
    samp_fit_args = list(fit_type = "external"),
    samp_fit_ext = ipcw_est_sl,
    g_exp_fit_args = list(fit_type = "external"),
    gn_exp_fit_ext = gn_exp_sl[C_samp == 1, ],
    g_cens_fit_args = list(fit_type = "external"),
    gn_cens_fit_ext = gn_cens_sl[C_samp == 1],
    Q_fit_args = list(fit_type = "external"),
    Qn_fit_ext = Qn_est_sl[C_samp == 1, ],
    eif_reg_type = "hal"
  )

  test_that("IPCW-TMLE and IPCW-one-step match with SL for sampling", {
    ipcw_tmle_sl_psi <- as.numeric(ipcw_tmle_sl$psi)
    ipcw_os_sl_psi <- as.numeric(ipcw_os_sl$psi)
    expect_equal(ipcw_tmle_sl_psi, ipcw_os_sl_psi, tol = 1e-2)
  })
}
