#' Efficient Influence Function of the Shifted Continuous Treatment Parameter
#'
#' The efficient influence function of the target parameter takes the form
#' D(P)(o) = H(a,w)(y - \bar{Q}(a,w)) + \bar{Q}(d(a,w)) - \psi(P)
#'
#' @param Y A \code{numeric} vector of observed outcomes.
#' @param A A \code{numeric} vector of observed treatments.
#' @param W A \code{matrix} or \code{data.frame} of baseline covariates.
#' @param Qn Function to compute the outcome regression: Q(A, W) = E(Y | A, W).
#' @param gn Function to compute the propensity score: g(A, W) = density(A | W).
#' @param delta A \code{numeric} for the shift to be placed on the treatment of
#' interest (i.e., the effect of shifting treatment \code{A} by delta units).
#' @param tol A \code{numeric} for the tolerance for measuring convergence of
#' the parametric fluctuation.
#' @param iter_max A \code{numeric} for the maximum number of iterations.
#' @param A_val A \code{vector} of \code{numeric} values for the points in the
#' range of the treatment \code{A} to approximate integrals by Riemmann sums.
#' Note that these must be equally spaced along a grid.
#'
#' @importFrom stats var
#'
#' @export
#'
#' @author Nima Hejazi
#' @author David Benkeser
#
eif <- function(Y, Hn, Qn_mat, fluct_fit) {
    # extract outcome regression for the 
    Qn_dAW <- Qn_mat$dAW
    Qn_AW <- Qn_mat$AW

    # maybe some other stuff might happen...

    # fit parameter estimate
    Qn_star <- predict(fluct_fit, newdata = data.frame(Qn_dAW),
                       type = "response")
    est_psi <- mean(Qn_star)

    # compute the efficient influence function
    eif <- Hn * (Y - Qn_AW) + Qn_dAW - est_psi

    # find the variance of the parameter
    var_psi <- mean(eif^2)

    # bundle EIC for the output
    out <- list(eif = eif, var = var_psi)
    return(out)
}

