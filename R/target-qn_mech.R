#' Targeting of Weighted Censoring Distribution
#'
#' DESCRIPTION HERE
#'
#' @param Qn_shift ...
#' @param ipc_weights ...
#' @param data_in ...
#'
#' @importFrom dplyr add_count select filter "%>%"
#' @importFrom tibble as_tibble
#' @importFrom rootSolve uniroot.all
#'
#' @keywords internal
#'
#' @export
#'
#' @author Nima Hejazi
#' @author David Benkeser
#
target_qn <- function(Qn_shift, ipc_weights, data_in) {
    # compute qn as the scaled sum of similar observations times IPC weights
    ind_a_w <- data_in %>%
      dplyr::add_count(W, round(A)) %>%
      dplyr::select(n) %>%
      unlist() %>%
      as.numeric()
    ipcw_trimmed <- ipc_weights %>%
      tibble::as_tibble() %>%
      dplyr::filter(value > 0)
    qn_dens <- ((ind_a_w * ipcw_trimmed) / nrow(data_in)) %>%
      unlist() %>%
      as.numeric()

    # compute D^F_W part of the EIF
    Dfw_eif <- Qn_shift - Qn_shift * qn_dens

    # bookkeeping with zeros for both D^F_W and qn terms
    qn_with_zeros <- rep(0, length(ipc_weights))
    Dfw_with_zeros <- rep(0, length(ipc_weights))
    qn_with_zeros[ipc_weights != 0] <- qn_dens
    Dfw_with_zeros[ipc_weights != 0] <- Dfw_eif

    # solve for root of the score equation of epsilon
    eps_score <- function(eps_in, Dfw, qn_term, ipc_weights) {
      eps_vec <- rep(eps_in, length(ipc_weights))
      mean(ipc_weights * ((Dfw * qn_term) / (1 + eps_vec * (Dfw * qn_term))))
    }
    eps_n <- rootSolve::uniroot.all(eps_score, interval = c(-0.5, 0.5),
                                    Dfw = Dfw_with_zeros,
                                    qn_term = qn_with_zeros,
                                    ipc_weights = ipc_weights)
    return(eps_n)
}

