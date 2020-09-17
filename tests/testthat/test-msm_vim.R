context("Marginal structural model summarizes TMLE and one-step")
library(data.table)
set.seed(172943)

if (require("sl3")) {
  # Example based on the data-generating mechanism presented in the simulation
  n <- 100
  W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
  A <- rpois(n, lambda = exp(3 + .3 * log(W$W1) - 0.2 * exp(W$W1) * W$W2))
  Y <- -1 + 0.05 * A - 0.02 * A * W$W2 + 0.2 * A * tan(W$W1^2) - 0.02 * W$W1 *
    W$W2 + 0.1 * A * W$W1 * W$W2 + rnorm(n, 0, 0.1)
  C <- rbinom(n, 1, plogis(rowSums(W) + Y))
  delta_grid <- seq(-0.5, 0.5, 0.5)
  EY <- mean(Y)
  EY_delta <- lapply(delta_grid, function(delta) {
    mean(-1 + 0.05 * A - 0.02 * (A + delta) * W$W2 + 0.2 * A * tan(W$W1^2) -
      0.02 * W$W1 * W$W2 + 0.1 * A * W$W1 * W$W2)
  })
  msm_lm <- lm.fit(y = do.call(c, EY_delta), x = cbind(1, delta_grid))

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
  sl_density <- Lrnr_sl$new(
    learners = g_lib,
    metalearner = Lrnr_solnp_density$new()
  )

  # fit MSM over TMLEs
  tmle_fits <- lapply(delta_grid, function(delta) {
    tmle <- txshift(
      Y = Y, A = A, W = W, delta = delta, estimator = "tmle",
      g_fit_args = list(fit_type = "sl", sl_learners_density = sl_density),
      Q_fit_args = list(fit_type = "sl", sl_learners = sl)
    )
    return(tmle)
  })
  tmle_psi <- do.call(c, lapply(tmle_fits, `[[`, "psi"))

  msm_tmle <- msm_vimshift(
    W = W, A = A, Y = Y,
    delta_grid = delta_grid, estimator = "tmle", ci_type = "marginal",
    g_fit_args = list(fit_type = "sl", sl_learners_density = sl_density),
    Q_fit_args = list(fit_type = "sl", sl_learners = sl)
  )

  test_that("TML point estimates match with MSM-based estimates", {
    expect_equal(tmle_psi, msm_tmle$param_est$psi, tol = 1e-2)
  })

  # fit MSM over one-step estimates
  os_fits <- lapply(delta_grid, function(delta) {
    os <- txshift(
      Y = Y, A = A, W = W, delta = delta, estimator = "onestep",
      g_fit_args = list(fit_type = "sl", sl_learners_density = sl_density),
      Q_fit_args = list(fit_type = "sl", sl_learners = sl)
    )
    return(os)
  })
  os_psi <- do.call(c, lapply(os_fits, `[[`, "psi"))

  msm_os <- msm_vimshift(
    W = W, A = A, Y = Y, estimator = "onestep", ci_type = "marginal",
    g_fit_args = list(fit_type = "sl", sl_learners_density = sl_density),
    Q_fit_args = list(fit_type = "sl", sl_learners = sl),
  )

  test_that("One-step point estimates match with MSM-based estimates", {
    expect_equal(os_psi, msm_os$param_est$psi, tol = 1e-2)
  })
}
