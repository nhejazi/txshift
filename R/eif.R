#' Efficient Influence Function of the Shifted Continuous Treatment Parameter
#'
#' The efficient influence function of the target parameter takes the form
#' D(P)(o) = H(a,w)(y - \bar{Q}(a,w)) + \bar{Q}(d(a,w)) - \psi(P)
#'
#' @param Y A \code{numeric} vector of observed outcomes.
#' @param Hn ...
#' @param Qn_mat ...
#' @param fluct_fit Model fit object...
#'
#' @importFrom stats var
#'
#' @keywords internal
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

