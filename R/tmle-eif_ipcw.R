#' Compute IPCW Part of the Efficient Influence Function for IPCW-TMLE
#'
#' Computes an additional component of the efficient influence function needed
#' for efficient estimation of IPCW-TMLEs. This takes the form:
#' %0=P_n(\Delta-\Pi_n^*(V)))\frac{\E(D^F(P^0_{X,n})\mid\Delta=1,V)}{\Pi^*_n(V)}
#'
#' @param fluc_fit_out ...
#' @param Hn ...
#' @param Y A \code{numeric} vector of observed outcomes.
#' @param ipc_weights ...
#' @param tol_eif ...
#'
#' @importFrom stats var
#' @importFrom dplyr if_else
#'
#' @keywords internal
#'
#' @export
#'
#' @author Nima Hejazi
#' @author David Benkeser
#
ipcwtmle_eif <- function(fluc_fit_out,
                         Hn,
                         Y,
                         C,
                         ipc_weights = rep(1, length(Y)),
                         tol_eif = 1e-7) {
  # compute TMLE
  psi <- mean(ipc_weights * fluc_fit_out$Qn_shift_star)

  # compute the efficient influence function (EIF) / canonical gradient (D*)
  eif <- ipc_weights * Hn$noshift * (Y - fluc_fit_out$Qn_noshift_star) +
    ipc_weights * (fluc_fit_out$Qn_shift_star - psi)


  # compute the variance based on the EIF
  # NOTE: scale by length of observations to get on same scale as parameter
  # NOTE: var(eif) and mean(eif^2) are nearly equivalent
  # var_eif <- mean(eif^2) / length(Y)
  var_eif <- stats::var(eif) / length(Y)

  # return the variance and the EIF value at each observation
  out <- list(psi = psi, var = var_eif, eif = eif)
  return(out)
}

