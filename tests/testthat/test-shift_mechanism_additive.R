context("additive treatment shifting mechanism works as expected")
library(data.table)
set.seed(73294)

# Example based on the data-generating mechanism presented in the simulation
n <- 100
W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
A <- rpois(n, lambda = exp(3 + .3 * log(W$W1) - 0.2 * exp(W$W1) * W$W2))
Y <- rbinom(
  n, 1,
  plogis(-1 + 0.05 * A - 0.02 * A * W$W2 + 0.2 * A * tan(W$W1^2) -
    0.02 * W$W1 * W$W2 + 0.1 * A * W$W1 * W$W2)
)
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

gn_spec_fitted <- as.data.table(
  lapply(
    c(-delta_shift, 0, delta_shift, 2 * delta_shift),
    function(delta) {
      txshift:::shift_additive(A = A, delta = delta)
    }
  )
)
setnames(gn_spec_fitted, c("downshift", "noshift", "upshift", "upupshift"))

test_that("Simple shifting function shifts downward by magnitude delta", {
  expect_equal(unique(gn_spec_fitted$downshift - A), -delta_shift)
})

test_that("Simple shifting function does not induce shift when delta = 0", {
  expect_equal(unique(gn_spec_fitted$noshift - A), 0)
})

test_that("Simple shifting function shifts upward by magnitude delta", {
  expect_equal(unique(gn_spec_fitted$upshift - A), delta_shift)
})

test_that("Simple shifting function shifts upward by twice magnitude delta", {
  expect_equal(unique(gn_spec_fitted$upupshift - A), 2 * delta_shift)
})
