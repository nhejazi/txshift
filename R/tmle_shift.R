#' Estimate the effect of a continuous treatment
#'
#' @param Y Outcome vector.
#' @param A Treatment vector.
#' @param W Covariate \code{matrix} or \code{data.frame}.
#' @param Qn Function to compute Q(A, W) = E(Y | A, W).
#' @param gn Function to compute g(A, W) = density(A | W)
#' @param delta Shift value of interest (i.e., compute the effect of shifting
#' the treatment A by delta units.
#' @param tol Convergence tolerance for parametric fluctuation.
#' @param iter_max Maximum number of iterations.
#' @param A_val Points in the range of the treatment A to approximate integrals
#' by Riemmann sums. These must be equally spaced along a grid.
#'
#' @importFrom stats var
#'
#' @export
#'
#' @author Ivan Diaz
#' @author Nima Hejazi
#'
#' @examples
#' n <- 100
#' W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
#' A <- rpois(n, lambda = exp(3 + 0.3 * log(W$W1) - 0.2 * exp(W$W1) * W$W2))
#' Y <- rbinom(n, 1, plogis(-1 + 0.05 * A - 0.02 * A * W$W2 +
#'                          0.2 * A * tan(W$W1^2) - 0.02 * W$W1 * W$W2 +
#'                          0.1 * A * W$W1 * W$W2))
#' fitA.0 <- glm(A ~ I(log(W1)) + I(exp(W1)):W2, family = poisson,
#'               data = data.frame(A, W))
#' fitY.0 <- glm(Y ~ A + A:W2 + A:I(tan(W1^2)) + W1:W2 + A:W1:W2,
#'               family = binomial, data = data.frame(A, W))
#' gn.0  <- function(A = A, W = W) {
#'   dpois(A, lambda = predict(fitA.0, newdata = W, type = "response"))
#' }
#' Qn.0 <- function(A = A, W = W) {
#'   predict(fitY.0, newdata = data.frame(A, W, row.names = NULL),
#'           type = "response")
#' }
#' tmle00 <- tmle_shift(Y = Y, A = A, W = W, Qn = Qn.0, gn = gn.0, delta = 2,
#'                      tol = 1e-4, iter_max = 5, A_val = seq(1, 60, 1))
#
tmle_shift <- function(Y, A, W, Qn, gn, delta, tol = 1e-5, iter_max = 5, A_val) {

  # interval partition length, A_val assumed equally spaced
  h_int <- A_val[3] - A_val[2]

  # this function takes as input initial estimator of Q and g and returns
  # their updated value
  ini.out <- f_iter(Qn = Qn, gn = gn, gn0d = NULL, prev_sum = 0, first = TRUE,
                    h_int = h_int, W = W, A = A, A_val = A_val, Y = Y,
                    delta = delta)

  gn0d <- ini.out$gn0d
  iter = 0

  # iterative procedure
  while(abs(ini.out$eps) > tol & iter <= iter_max){
    iter = iter + 1
    new.out <- f_iter(Qn = ini.out$Qn, gn = ini.out$gn, gn0d = gn0d,
                      prev_sum = ini.out$prev_sum, first = FALSE, h_int = h_int,
                      W = W, A = A, A_val = A_val, Y = Y, delta = delta)
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

#' Fluctuation (iterative)
#'
#' @param Qn Function to compute Q(A, W) = E(Y | A, W).
#' @param gn Function to compute g(A, W) = density(A | W).
#' @param gn0d \code{Numeric} computed for gn from the first iteration only.
#' @param prev_sum ...
#' @param first ...
#' @param W ...
#' @param A_val ...
#' @param h_int ...
#' @param A ...
#' @param Y ...
#' @param delta ...
#'
#' @importFrom stats uniroot
#'
#' @author Ivan Diaz
#' @author Nima Hejazi
#
f_iter <- function(Qn, gn, gn0d = NULL, prev_sum = 0, first = FALSE, h_int,
                   W, A, A_val, Y, delta) {

  # numerical integrals and equation (7)
  Qnd <- t(sapply(seq_len(nrow(W)), function(i) Qn(A_val + delta, W[i, ])))
  gnd <- t(sapply(seq_len(nrow(W)), function(i) gn(A_val, W[i, ])))
  gnd <- gnd / rowSums(gnd)
  if (first) gn0d <- gnd

  EQnd <- rowSums(Qnd * gnd) * h_int
  D2   <- Qnd - EQnd
  QnAW <- Qn(A, W)
  H1   <- gn(A - delta, W) / gn(A, W)

  # equation (8)
  est_eqn_min  <- stats::uniroot(est_eqn, c(-1, 1),  Y = Y, A = A, W = W,
                                 delta = delta, QnAW = QnAW, Qn = Qn, H1 = H1,
                                 gn0d = gn0d, EQnd = EQnd, D2 = D2,
                                 prev_sum = prev_sum)
  eps <- est_eqn_min$root

  # updated values
  gn.new   <- function(a, w) exp(eps * Qn(a + delta, w)) * gn(a, w)
  Qn.new   <- function(a, w) Qn(a, w) + eps * gn(a - delta, w) / gn(a, w)
  prev_sum <- prev_sum + eps * D2
  return(list(Qn = Qn.new, gn = gn.new, prev_sum = prev_sum, eps = eps,
              gn0d = gn0d))
}

################################################################################

#' Estimating equation
#'
#' @param eps ...
#' @param Y ...
#' @param A ...
#' @param W ...
#' @param delta ...
#' @param QnAW ...
#' @param Qn ...
#' @param H1 ...
#' @param gn0d ...
#' @param prev_sum ...
#' @param EQnd ...
#' @param D2 ...
#'
#' @author Ivan Diaz
#' @author Nima Hejazi
#
est_eqn <- function(eps, Y, A, W, delta, QnAW, Qn, H1, gn0d, EQnd, D2,
                    prev_sum) {

  sum((Y - (QnAW + eps * H1)) * H1 + (Qn(A + delta, W) - EQnd) -
      rowSums(D2 * exp(eps * D2 + prev_sum) * gn0d) /
      rowSums(exp(eps * D2 + prev_sum) * gn0d))
}

