context("One-step and TML estimator implementations agree")
library(data.table)
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
gn_spec_fitted <- as.data.table(
  lapply(
    c(-delta_shift, 0, delta_shift, 2 * delta_shift),
    function(shift_value) {
      gn_out <- gn.0(A = A + shift_value, W = W)
    }
  )
)
setnames(gn_spec_fitted, c("downshift", "noshift", "upshift", "upupshift"))

# NOTE: should also use true Q for good measure (truth includes interactions)
Qn_spec_fitted <- as.data.table(
  lapply(c(0, delta_shift), function(shift_value) {
    Qn_out <- Qn.0(A = A + shift_value, W = W)
  })
)
setnames(Qn_spec_fitted, c("noshift", "upshift"))

# fit TMLE
tmle_txshift <- txshift(
  Y = Y, A = A, W = W, delta = delta_shift,
  g_fit = list(fit_type = "fit_spec"),
  Q_fit = list(fit_type = "fit_spec"),
  gn_fit_spec = gn_spec_fitted,
  Qn_fit_spec = Qn_spec_fitted,
  estimator = "tmle"
)
tmle_psi <- as.numeric(tmle_txshift$psi)

# fit one-step
os_txshift <- txshift(
  Y = Y, A = A, W = W, delta = delta_shift,
  g_fit = list(fit_type = "fit_spec"),
  Q_fit = list(fit_type = "fit_spec"),
  gn_fit_spec = gn_spec_fitted,
  Qn_fit_spec = Qn_spec_fitted,
  estimator = "onestep"
)
os_psi <- as.numeric(os_txshift$psi)

# test for reasonable equality between estimators
test_that("TMLE and one-step implementations match closely", {
  expect_equal(tmle_psi, os_psi, tol = 1e-3)
})

# Example of IPCW-based estimators by adding censoring node
ipcw_tmle <- txshift(
  W = W, A = A, Y = Y, delta = delta_shift,
  C = C, V = c("W", "Y"),
  estimator = "tmle",
  max_iter = 5,
  ipcw_fit_args = list(fit_type = "glm"),
  g_fit = list(fit_type = "fit_spec"),
  Q_fit = list(fit_type = "fit_spec"),
  gn_fit_spec = gn_spec_fitted[C == 1, ],
  Qn_fit_spec = Qn_spec_fitted[C == 1, ],
  eif_reg_type = "glm"
)
ipcw_tmle_psi <- as.numeric(ipcw_tmle$psi)

ipcw_os <- txshift(
  W = W, A = A, Y = Y, delta = delta_shift,
  C = C, V = c("W", "Y"),
  estimator = "onestep",
  ipcw_fit_args = list(fit_type = "glm"),
  g_fit = list(fit_type = "fit_spec"),
  Q_fit = list(fit_type = "fit_spec"),
  gn_fit_spec = gn_spec_fitted[C == 1, ],
  Qn_fit_spec = Qn_spec_fitted[C == 1, ],
  eif_reg_type = "glm"
)
ipcw_os_psi <- as.numeric(ipcw_os$psi)

# test for reasonable equality between estimators
test_that("TMLE and one-step match reasonably well under censoring", {
  expect_equal(ipcw_tmle_psi, ipcw_os_psi, tol = 1e-3)
})
