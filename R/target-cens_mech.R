#' Targeting of Censoring Mechanism Estimate
#'
#' DESCRIPTION HERE
#'
#' @param eps_n ...
#' @param Q_n ...
#' @param Pi_n ...
#' @param delta_ind ...
#'
#' @keywords internal
#'
#' @export
#'
#' @author Nima Hejazi
#' @author David Benkeser
#
target_pi <- function(eps_n, Q_n, Pi_n, delta_ind) {
    # first, define a part of the targeting equation that appears twice
    Qn_diff = Q_n - mean((delta_ind / Pi_n) * Q_n)
    eps_score <- mean((delta_ind / Pi_n) * (Qn_diff / (1 + eps_n(Qn_diff))))
    c(eps_new = eps_score)
}

