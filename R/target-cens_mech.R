#' Targeting of Censoring Mechanism Estimate
#'
#' DESCRIPTION HERE
#'
#' @param Qn ...
#' @param ipc_weights ...
#'
#' @importFrom rootSolve uniroot.all
#'
#' @keywords internal
#'
#' @export
#'
#' @author Nima Hejazi
#' @author David Benkeser
#
target_pi <- function(Qn, ipc_weights) {
    # bookkeeping with zeros in the censoring mechanism
    Qn_with_zeros <- rep(0, length(ipc_weights))
    Qn_with_zeros[ipc_weights != 0] <- Qn

    # compute IPC-weighted outcome while preserving zeros
    Qn_weighted <- ipc_weights * Qn_with_zeros

    # first, define a part of the targeting equation that appears twice
    Qn_weighted_updated <- Qn_with_zeros - mean(Qn_weighted)

    # solve for root of the score equation of epsilon
    eps_score <- function(eps_in, Qn_with_weights, ipc_weights) {
      eps_vec <- rep(eps_in, length(ipc_weights))
      mean(ipc_weights * (Qn_with_weights / (1 + eps_vec * Qn_with_weights)))
    }
    eps_n <- rootSolve::uniroot.all(eps_score, interval = c(-1, 1),
                                    Qn_with_weights = Qn_weighted_updated,
                                    ipc_weights = ipc_weights)
    return(eps_n)
}

