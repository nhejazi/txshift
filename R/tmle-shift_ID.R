#' TML Estimate of the Effect of a Continuous Treatment
#'
#' @param Y A \code{numeric} vector of observed outcomes.
#' @param A A \code{numeric} vector of observed treatments.
#' @param W A \code{matrix} or \code{data.frame} of baseline covariates.
#' @param Qn Function to compute the outcome regression: Q(A, W) = E(Y | A, W).
#' @param gn Function to compute the propensity score: g(A, W) = density(A | W).
#' @param delta A \code{numeric} for the shift to be placed on the treatment of
#'  interest (i.e., the effect of shifting treatment \code{A} by delta units).
#' @param tol A \code{numeric} for the tolerance for measuring convergence of
#'  the parametric fluctuation.
#' @param iter_max A \code{numeric} for the maximum number of iterations.
#' @param A_val A \code{vector} of \code{numeric} values for the points in the
#'  range of the treatment \code{A} to approximate integrals by Riemmann sums.
#'  Note that these must be equally spaced along a grid.
#'
#' @importFrom stats var
#'
#' @keywords internal
#'
#' @export
#'
#' @author Iván Díaz
#' @author Nima Hejazi
#'
#' @examples
#' n <- 100
#' W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
#' A <- rpois(n, lambda = exp(3 + 0.3 * log(W$W1) - 0.2 * exp(W$W1) * W$W2))
#' Y <- rbinom(n, 1, plogis(-1 + 0.05 * A - 0.02 * A * W$W2 +
#'   0.2 * A * tan(W$W1^2) - 0.02 * W$W1 * W$W2 +
#'   0.1 * A * W$W1 * W$W2))
#' fitA.0 <- glm(A ~ I(log(W1)) + I(exp(W1)):W2,
#'   family = poisson,
#'   data = data.frame(A, W)
#' )
#' fitY.0 <- glm(Y ~ A + A:W2 + A:I(tan(W1^2)) + W1:W2 + A:W1:W2,
#'   family = binomial, data = data.frame(A, W)
#' )
#' gn.0 <- function(A = A, W = W) {
#'   dpois(A, lambda = predict(fitA.0, newdata = W, type = "response"))
#' }
#' Qn.0 <- function(A = A, W = W) {
#'   predict(fitY.0,
#'     newdata = data.frame(A, W, row.names = NULL),
#'     type = "response"
#'   )
#' }
#' tmle00 <- tmle_shift_orig(
#'   Y = Y, A = A, W = W, Qn = Qn.0, gn = gn.0, delta = 2,
#'   tol = 1e-4, iter_max = 5, A_val = seq(1, 60, 1)
#' )
#' #
tmle_shift_orig <- function(Y, A, W,
                            Qn, gn,
                            delta, tol = 1e-5, iter_max = 5, A_val) {

  # interval partition length, A_val assumed equally spaced
  h_int <- A_val[3] - A_val[2]

  # this function takes as input initial estimator of Q and g and returns
  # their updated value
  ini.out <- f_iter(
    Qn = Qn, gn = gn, gn0d = NULL, prev_sum = 0, first = TRUE,
    h_int = h_int, W = W, A = A, A_val = A_val, Y = Y,
    delta = delta
  )

  gn0d <- ini.out$gn0d
  iter <- 0

  # iterative procedure
  while (abs(ini.out$eps) > tol & iter <= iter_max) {
    iter <- iter + 1
    new.out <- f_iter(
      Qn = ini.out$Qn, gn = ini.out$gn, gn0d = gn0d,
      prev_sum = ini.out$prev_sum, first = FALSE, h_int = h_int,
      W = W, A = A, A_val = A_val, Y = Y, delta = delta
    )
    ini.out <- new.out
  }
  Qnd <- t(sapply(seq_len(nrow(W)), function(i) ini.out$Qn(A_val + delta, W[i, ])))
  gnd <- t(sapply(seq_len(nrow(W)), function(i) ini.out$gn(A_val, W[i, ])))
  gnd <- gnd / rowSums(gnd)

  # plug in tmle
  psi.hat <- mean(rowSums(Qnd * gnd) * h_int)

  # influence curve of tmle
  IC <- (Y - ini.out$Qn(A, W)) * ini.out$gn(A - delta, W) / ini.out$gn(A, W) +
    ini.out$Qn(A + delta, W) - psi.hat
  var.hat <- stats::var(IC) / length(Y)
  return(c(psi.hat = psi.hat, var.hat = var.hat, IC = IC))
}

################################################################################

#' Iterative Fluctuation
#'
#' @param Qn Function to compute the outcome regression: Q(A, W) = E(Y | A, W).
#' @param gn Function to compute the propensity score: g(A, W) = density(A | W).
#' @param gn0d A \code{Numeric} computed for the propensity score (gn) from the
#'  first iteration only.
#' @param prev_sum The sum computed from this fluctuation step in the previous
#'  iteration.
#' @param first A \code{logical} for whether this is the first iteration of the
#'  iterative fluctutation step.
#' @param h_int ???
#' @param Y A \code{numeric} vector of observed outcomes.
#' @param A A \code{numeric} vector of observed treatments.
#' @param W A \code{matrix} or \code{data.frame} of baseline covariates.
#' @param delta A \code{numeric} for the shift to be placed on the treatment of
#'  interest (i.e., the effect of shifting treatment \code{A} by delta units).
#' @param A_val A \code{vector} of \code{numeric} values for the points in the
#'  range of the treatment \code{A} to approximate integrals by Riemmann sums.
#'  Note that these must be equally spaced along a grid.
#'
#' @importFrom stats uniroot
#'
#' @author Iván Díaz
#' @author Nima Hejazi
#'
#' @keywords internal
#
f_iter <- function(Qn, gn, gn0d = NULL, prev_sum = 0, first = FALSE, h_int,
                   Y, A, W, delta, A_val) {

  # numerical integrals and equation (7)
  Qnd <- t(sapply(seq_len(nrow(W)), function(i) Qn(A_val + delta, W[i, ])))
  gnd <- t(sapply(seq_len(nrow(W)), function(i) gn(A_val, W[i, ])))
  gnd <- gnd / rowSums(gnd)
  if (first) gn0d <- gnd

  EQnd <- rowSums(Qnd * gnd) * h_int
  D2 <- Qnd - EQnd
  QnAW <- Qn(A, W)
  H1 <- gn(A - delta, W) / gn(A, W)

  # equation (8)
  est_eqn_min <- stats::uniroot(
    est_eqn, c(-1, 1),
    Y = Y, A = A, W = W,
    delta = delta, QnAW = QnAW, Qn = Qn, H1 = H1,
    gn0d = gn0d, EQnd = EQnd, D2 = D2,
    prev_sum = prev_sum
  )
  eps <- est_eqn_min$root

  # updated values
  gn_new <- function(a, w) exp(eps * Qn(a + delta, w)) * gn(a, w)
  Qn_new <- function(a, w) Qn(a, w) + eps * gn(a - delta, w) / gn(a, w)
  prev_sum <- prev_sum + eps * D2
  return(list(
    Qn = Qn_new, gn = gn_new, prev_sum = prev_sum, eps = eps,
    gn0d = gn0d
  ))
}

################################################################################

#' Estimating Equation
#'
#' @param eps ...
#' @param Qn Function to compute the outcome regression: Q(A, W) = E(Y | A, W).
#' @param QnAW ...
#' @param EQnd ...
#' @param gn0d A \code{Numeric} computed for the propensity score (gn) from the
#'  first iteration only.
#' @param H1 Ratio obtained from comparing the propensity score (gn) before and
#'  after application of the shift \code{delta}.
#' @param D2 ...
#' @param prev_sum The sum computed from this fluctuation step in the previous
#'  iteration.
#' @param Y A \code{numeric} vector of observed outcomes.
#' @param A A \code{numeric} vector of observed treatments.
#' @param W A \code{matrix} or \code{data.frame} of baseline covariates.
#' @param delta A \code{numeric} for the shift to be placed on the treatment of
#'  interest (i.e., the effect of shifting treatment \code{A} by delta units).
#'
#' @author Iván Díaz
#' @author Nima Hejazi
#'
#' @keywords internal
#
est_eqn <- function(eps, QnAW, Qn, H1, gn0d, EQnd, D2, prev_sum, Y, A, W,
                    delta) {
  sum((Y - (QnAW + eps * H1)) * H1 + (Qn(A + delta, W) - EQnd) -
    rowSums(D2 * exp(eps * D2 + prev_sum) * gn0d) /
      rowSums(exp(eps * D2 + prev_sum) * gn0d))
}
