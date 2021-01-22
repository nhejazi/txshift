context("Implementations of TML estimator variants agree")
library(data.table)
library(rlang)
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
delta_shift <- 2

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

# NOTE: using true density like Ivan does
gn_ext_fitted <- as.data.table(
  lapply(
    c(-delta_shift, 0, delta_shift, 2 * delta_shift),
    function(shift_value) {
      gn_out <- gn.0(A = A + shift_value, W = W)
    }
  )
) %>% set_names(c("downshift", "noshift", "upshift", "upupshift"))

# NOTE: should also use true Q for good measure (truth includes interactions)
Qn_ext_fitted <- as.data.table(
  lapply(c(0, delta_shift), function(shift_value) {
    Qn_out <- Qn.0(A = A + shift_value, W = W)
  })
) %>% set_names(c("noshift", "upshift"))

# fit TMLE
tmle_txshift_std <- txshift(
  Y = Y, A = A, W = W, delta = delta_shift,
  g_exp_fit_args = list(fit_type = "external"),
  gn_exp_fit_ext = gn_ext_fitted,
  Q_fit_args = list(fit_type = "external"),
  Qn_fit_ext = Qn_ext_fitted,
  estimator = "tmle",
  fluctuation = "standard"
)
psi_std <- as.numeric(tmle_txshift_std$psi)

# fit one-step
tmle_txshift_wts <- txshift(
  Y = Y, A = A, W = W, delta = delta_shift,
  g_exp_fit_args = list(fit_type = "external"),
  gn_exp_fit_ext = gn_ext_fitted,
  Q_fit_args = list(fit_type = "external"),
  Qn_fit_ext = Qn_ext_fitted,
  estimator = "tmle",
  fluctuation = "weighted"
)
psi_wts <- as.numeric(tmle_txshift_wts$psi)

# test for reasonable equality between estimators
test_that("TMLEs match for weighted v. standard fluctuation", {
  expect_equal(psi_std, psi_wts, tol = 1e-3)
})

# Example of IPCW-based estimators by adding censoring node
ipcw_tmle_glm_std <- txshift(
  W = W, A = A, Y = Y, delta = delta_shift,
  C_samp = C, V = c("W", "Y"),
  estimator = "tmle",
  max_iter = 5,
  samp_fit_args = list(fit_type = "glm"),
  g_exp_fit_args = list(fit_type = "external"),
  gn_exp_fit_ext = gn_ext_fitted[C == 1, ],
  Q_fit_args = list(fit_type = "external"),
  Qn_fit_ext = Qn_ext_fitted[C == 1, ],
  eif_reg_type = "glm",
  fluctuation = "standard"
)
psi_ipcw_tmle_glm_std <- as.numeric(ipcw_tmle_glm_std$psi)

ipcw_tmle_glm_wts <- txshift(
  W = W, A = A, Y = Y, delta = delta_shift,
  C_samp = C, V = c("W", "Y"),
  estimator = "tmle",
  max_iter = 5,
  samp_fit_args = list(fit_type = "glm"),
  g_exp_fit_args = list(fit_type = "external"),
  gn_exp_fit_ext = gn_ext_fitted[C == 1, ],
  Q_fit_args = list(fit_type = "external"),
  Qn_fit_ext = Qn_ext_fitted[C == 1, ],
  eif_reg_type = "glm",
  fluctuation = "weighted"
)
psi_ipcw_tmle_glm_wts <- as.numeric(ipcw_tmle_glm_wts$psi)

# test for reasonable equality between estimators
test_that("IPCW-TMLEs w/ EIF-GLM match for weighted v. standard fluctuation", {
  expect_equal(psi_ipcw_tmle_glm_std, psi_ipcw_tmle_glm_wts, tol = 1e-3)
})

# Example of IPCW-based estimators by adding censoring node
ipcw_tmle_hal_std <- txshift(
  W = W, A = A, Y = Y, delta = delta_shift,
  C_samp = C, V = c("W", "Y"),
  estimator = "tmle",
  max_iter = 5,
  samp_fit_args = list(fit_type = "glm"),
  g_exp_fit_args = list(fit_type = "external"),
  gn_exp_fit_ext = gn_ext_fitted[C == 1, ],
  Q_fit_args = list(fit_type = "external"),
  Qn_fit_ext = Qn_ext_fitted[C == 1, ],
  eif_reg_type = "hal",
  fluctuation = "standard"
)
psi_ipcw_tmle_hal_std <- as.numeric(ipcw_tmle_hal_std$psi)

ipcw_tmle_hal_wts <- txshift(
  W = W, A = A, Y = Y, delta = delta_shift,
  C_samp = C, V = c("W", "Y"),
  estimator = "tmle",
  max_iter = 5,
  samp_fit_args = list(fit_type = "glm"),
  g_exp_fit_args = list(fit_type = "external"),
  gn_exp_fit_ext = gn_ext_fitted[C == 1, ],
  Q_fit_args = list(fit_type = "external"),
  Qn_fit_ext = Qn_ext_fitted[C == 1, ],
  eif_reg_type = "hal",
  fluctuation = "weighted"
)
psi_ipcw_tmle_hal_wts <- as.numeric(ipcw_tmle_hal_wts$psi)

# test for reasonable equality between estimators
test_that("IPCW-TMLEs w/ EIF-HAL match for weighted v. standard fluctuation", {
  expect_equal(psi_ipcw_tmle_hal_std, psi_ipcw_tmle_hal_wts, tol = 1e-3)
})
