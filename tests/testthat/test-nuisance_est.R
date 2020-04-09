context("Estimation of nuisance parameters")
library(data.table)
library(sl3)
set.seed(172943)

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
V <- as.data.table(list(W, Y = Y))
delta_shift <- 2

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
sl_density <- Lrnr_sl$new(learners = g_lib,
                          metalearner = Lrnr_solnp_density$new())

# fit exposure mechanism with HAL, SL, or GLM
gn_est_hal <- est_g(A = A, W = W,
                    delta = delta_shift,
                    fit_type = "hal",
                    haldensify_args = list(
                      n_bins = 5,
                      grid_type = "equal_mass",
                      lambda_seq = exp(seq(-1, -13, length = 300)),
                      use_future = FALSE)
                   )
gn_est_sl <- est_g(A = A, W = W,
                   delta = delta_shift,
                   fit_type = "sl",
                   sl_learners_density = sl_density)
gn_est_glm <- est_g(A = A, W = W,
                   delta = delta_shift,
                   fit_type = "sl",
                   sl_learners_density = hse_learner)

# fit outcome mechanism with GLM or SL
Qn_est_glm <- est_Q(Y = Y, A = A, W = W, delta = delta_shift,
                    fit_type = "glm", glm_formula = "Y ~ .")
Qn_est_sl <- est_Q(Y = Y, A = A, W = W, delta = delta_shift,
                   fit_type = "sl", sl_learners = sl)

# fit two-phase censoring mechanism with GLM or SL
ipcw_est_glm <- est_ipcw(V = V, Delta = C, fit_type = "glm")
ipcw_est_sl <- est_ipcw(V = V, Delta = C, fit_type = "sl", sl_learners = sl)

# test for errors when arguments are set inconsistently
test_that("SL-based nuisance estimation fails without SL library.", {
  expect_error(est_g(A = A, W = W, delta = delta_shift, fit_type = "sl"))
  expect_error(est_Q(A = A, W = W, delta = delta_shift, fit_type = "sl"))
  expect_error(est_ipcw(A = A, W = W, delta = delta_shift, fit_type = "sl"))
})

# fit TMLE and one-step for HAL and SL
tmle_ml <- txshift(
  Y = Y, A = A, W = W, delta = delta_shift,
  g_fit = list(fit_type = "external"),
  gn_fit_ext = gn_est_hal,
  Q_fit = list(fit_type = "external"),
  Qn_fit_ext = Qn_est_sl,
  estimator = "tmle"
)
os_ml <- txshift(
  Y = Y, A = A, W = W, delta = delta_shift,
  g_fit = list(fit_type = "external"),
  gn_fit_ext = gn_est_hal,
  Q_fit = list(fit_type = "external"),
  Qn_fit_ext = Qn_est_sl,
  estimator = "onestep"
)

test_that("TMLE and one-step match for HAL/SL nuisance configurations", {
  tmle_ml_psi <- as.numeric(tmle_ml$psi)
  os_ml_psi <- as.numeric(os_ml$psi)
  expect_equal(tmle_ml_psi, os_ml_psi, tol = 1e-2)
})

# fit TMLE and one-step for GLMs
tmle_glm <- txshift(
  Y = Y, A = A, W = W, delta = delta_shift,
  g_fit = list(fit_type = "external"),
  gn_fit_ext = gn_est_glm,
  Q_fit = list(fit_type = "external"),
  Qn_fit_ext = Qn_est_glm,
  estimator = "tmle"
)
os_glm <- txshift(
  Y = Y, A = A, W = W, delta = delta_shift,
  g_fit = list(fit_type = "external"),
  gn_fit_ext = gn_est_glm,
  Q_fit = list(fit_type = "external"),
  Qn_fit_ext = Qn_est_glm,
  estimator = "onestep"
)

test_that("TMLE and one-step match for GLM nuisance configurations", {
  tmle_glm_psi <- as.numeric(tmle_glm$psi)
  os_glm_psi <- as.numeric(os_glm$psi)
  expect_equal(tmle_glm_psi, os_glm_psi, tol = 1e-2)
})

# fit IPCW-TMLE and IPCW-one-step with SL for censoring estimation
ipcw_tmle_sl <- txshift(
  W = W, A = A, Y = Y, delta = delta_shift,
  C = C, V = c("W", "Y"),
  estimator = "tmle",
  max_iter = 5,
  ipcw_fit_args = list(fit_type = "external"),
  ipcw_fit_ext = ipcw_est_sl,
  g_fit = list(fit_type = "external"),
  gn_fit_ext = gn_est_hal[C == 1, ],
  Q_fit = list(fit_type = "external"),
  Qn_fit_ext = Qn_est_sl[C == 1, ],
  eif_reg_type = "hal"
)
ipcw_os_sl <- txshift(
  W = W, A = A, Y = Y, delta = delta_shift,
  C = C, V = c("W", "Y"),
  estimator = "onestep",
  ipcw_fit_args = list(fit_type = "external"),
  ipcw_fit_ext = ipcw_est_sl,
  g_fit = list(fit_type = "external"),
  gn_fit_ext = gn_est_hal[C == 1, ],
  Q_fit = list(fit_type = "external"),
  Qn_fit_ext = Qn_est_sl[C == 1, ],
  eif_reg_type = "hal"
)

test_that("IPCW-TMLE and IPCW-one-step match with SL censoring estimation", {
  ipcw_tmle_sl_psi <- as.numeric(ipcw_tmle_sl$psi)
  ipcw_os_sl_psi <- as.numeric(ipcw_os_sl$psi)
  expect_equal(ipcw_tmle_sl_psi, ipcw_os_sl_psi, tol = 1e-2)
})

# fit IPCW-TMLE and IPCW-one-step with GLM for censoring estimation
ipcw_tmle_glm <- txshift(
  W = W, A = A, Y = Y, delta = delta_shift,
  C = C, V = c("W", "Y"),
  estimator = "tmle",
  max_iter = 5,
  ipcw_fit_args = list(fit_type = "external"),
  ipcw_fit_ext = ipcw_est_glm,
  g_fit = list(fit_type = "external"),
  gn_fit_ext = gn_est_hal[C == 1, ],
  Q_fit = list(fit_type = "external"),
  Qn_fit_ext = Qn_est_sl[C == 1, ],
  eif_reg_type = "hal"
)
ipcw_os_glm <- txshift(
  W = W, A = A, Y = Y, delta = delta_shift,
  C = C, V = c("W", "Y"),
  estimator = "onestep",
  ipcw_fit_args = list(fit_type = "external"),
  ipcw_fit_ext = ipcw_est_glm,
  g_fit = list(fit_type = "external"),
  gn_fit_ext = gn_est_hal[C == 1, ],
  Q_fit = list(fit_type = "external"),
  Qn_fit_ext = Qn_est_sl[C == 1, ],
  eif_reg_type = "hal"
)

test_that("IPCW-TMLE and IPCW-one-step match with GLM censoring estimation", {
  ipcw_tmle_glm_psi <- as.numeric(ipcw_tmle_glm$psi)
  ipcw_os_glm_psi <- as.numeric(ipcw_os_glm$psi)
  expect_equal(ipcw_tmle_glm_psi, ipcw_os_glm_psi, tol = 1e-2)
})
