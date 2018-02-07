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
#' @importFrom stats uniroot
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
      dplyr::add_count(W, A) %>%
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
    Dfw_eif <- Qn_shift - sum(Qn_shift * qn_dens)

    # bookkeeping with zeros for both D^F_W and qn terms
    qn_with_zeros <- rep(0, length(ipc_weights))
    Dfw_with_zeros <- rep(0, length(ipc_weights))
    qn_with_zeros[ipc_weights != 0] <- qn_dens
    Dfw_with_zeros[ipc_weights != 0] <- Dfw_eif

    # solve for root of the score equation of epsilon
    eps_score <- function(eps_in, Dfw, ipc_weights) {
      out <- mean(ipc_weights * (Dfw  / (1 + eps_in * Dfw )))
    }
    # TODO: check that solution of epsilon is not on the boundary
    eps_n <- stats::uniroot(eps_score, interval = c(-1000, 1000),
                            Dfw = Dfw_with_zeros,
                            ipc_weights = ipc_weights)

    # update the estimated density qn based on epsilon
    qn_update <- (1 + eps_n$root * Dfw_with_zeros) * qn_with_zeros
    return(list(eps = eps_n$root, qn = qn_update))
}

