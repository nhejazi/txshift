context("tmle_shift_orig agrees with Diaz and van der Laan (2012)")
library(data.table)
library(rlang)
set.seed(73294)

################################################################################
## Original function from Diaz and van der Laan (2012), Biometrics
## https://github.com/idiazst/continuoustreat/blob/master/R/functions.R
################################################################################
tmle.shift <- function(Y, A, W, Qn, gn, delta, tol = 1e-5, iter.max = 5, Aval) {
  # interval partition length, Aval assumed equally spaced
  h.int <- Aval[3] - Aval[2]
  # this function takes as input initial estimator of Q and g and returns
  # their updated value
  f.iter <- function(Qn, gn, gn0d = NULL, prev.sum = 0, first = FALSE) {
    # numerical integrals and equation (7)
    Qnd <- t(sapply(1:nrow(W), function(i) Qn(Aval + delta, W[i, ])))
    gnd <- t(sapply(1:nrow(W), function(i) gn(Aval, W[i, ])))
    gnd <- gnd / rowSums(gnd)
    if (first) gn0d <- gnd
    EQnd <- rowSums(Qnd * gnd) * h.int
    D2 <- Qnd - EQnd
    QnAW <- Qn(A, W)
    H1 <- gn(A - delta, W) / gn(A, W)
    # equation (8)
    est.equation <- function(eps) {
      sum((Y - (QnAW + eps * H1)) * H1 + (Qn(A + delta, W) - EQnd) -
        rowSums(D2 * exp(eps * D2 + prev.sum) * gn0d) /
          rowSums(exp(eps * D2 + prev.sum) * gn0d))
    }
    eps <- uniroot(est.equation, c(-1, 1))$root
    # updated values
    gn.new <- function(a, w) exp(eps * Qn(a + delta, w)) * gn(a, w)
    Qn.new <- function(a, w) Qn(a, w) + eps * gn(a - delta, w) / gn(a, w)
    prev.sum <- prev.sum + eps * D2
    return(list(
      Qn = Qn.new, gn = gn.new, prev.sum =
        prev.sum, eps = eps, gn0d = gn0d
    ))
  }
  ini.out <- f.iter(Qn, gn, first = TRUE)
  gn0d <- ini.out$gn0d
  iter <- 0
  # iterative procedure
  while (abs(ini.out$eps) > tol & iter <= iter.max) {
    iter <- iter + 1
    new.out <- f.iter(ini.out$Qn, ini.out$gn, gn0d, ini.out$prev.sum)
    ini.out <- new.out
  }
  Qnd <- t(sapply(1:nrow(W), function(i) ini.out$Qn(Aval + delta, W[i, ])))
  gnd <- t(sapply(1:nrow(W), function(i) ini.out$gn(Aval, W[i, ])))
  gnd <- gnd / rowSums(gnd)
  # plug in tmle
  psi.hat <- mean(rowSums(Qnd * gnd) * h.int)
  # influence curve of tmle
  IC <- (Y - ini.out$Qn(A, W)) * ini.out$gn(A - delta, W) / ini.out$gn(A, W) +
    ini.out$Qn(A + delta, W) - psi.hat
  var.hat <- var(IC) / length(Y)
  return(c(psi.hat = psi.hat, var.hat = var.hat, IC = IC))
}

################################################################################

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

# run the two TMLE-shift algorithms
tmle_shift_2012 <- tmle.shift(
  Y = Y, A = A, W = W, Qn = Qn.0, gn = gn.0,
  delta = delta_shift, tol = 1e-4, iter.max = 5,
  Aval = seq(1, 60, 1)
)
tmle_2012_psi <- as.numeric(tmle_shift_2012[1])

# run the new txshift implementation of TMLE
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
tmle_txshift <- txshift(
  Y = Y, A = A, W = W, delta = delta_shift,
  g_fit = list(fit_type = "external"),
  gn_fit_ext = gn_ext_fitted,
  Q_fit = list(fit_type = "external"),
  Qn_fit_ext = Qn_ext_fitted
)
txshift_psi <- as.numeric(tmle_txshift$psi)

# test for equality between Ivan's modified code and txshift implementation
test_that("txshift implementation matches revised 2012 procedure closely", {
  expect_equal(tmle_2012_psi, txshift_psi, tol = 1e-3)
})
