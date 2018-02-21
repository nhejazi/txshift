#' Compute IPC-Weighted Component of Efficient Influence Function
#'
#' Computes an additional component of the efficient influence function needed
#' for efficient estimation of IPCW-TMLEs. This takes the form:
#' %0=P_n(\Delta-\Pi_n^*(V)))\frac{\E(D^F(P^0_{X,n})\mid\Delta=1,V)}{\Pi^*_n(V)}
#'
#' @param fluc_fit_out ...
#' @param data_in ...
#' @param Hn ...
#' @param Y A \code{numeric} vector of observed outcomes.
#' @param Delta ...
#' @param ipc_weights ...
#' @param tol_eif ...
#'
#' @importFrom stats var
#' @importFrom dplyr if_else add_count select mutate
#'
#' @keywords internal
#'
#' @export
#'
#' @author Nima Hejazi
#' @author David Benkeser
#
tmle_eif_ipcw <- function(fluc_fit_out,
                          data_in,
                          Hn,
                          Y,
                          Delta,
                          ipc_weights = rep(1, length(Y)),
                          tol_eif = 1e-7) {

  # compute TMLE of the treatment shift parameter
  param_obs_est <- rep(0, length(Delta))
  param_obs_est[Delta == 1] <- ipc_weights * fluc_fit_out$Qn_shift_star
  psi <- mean(param_obs_est)

  # compute the efficient influence function (EIF) / canonical gradient (D*)
  eif <- rep(0, length(Delta))
  eif[Delta == 1] <- ipc_weights * Hn$noshift *
    (Y - fluc_fit_out$Qn_noshift_star) +
    ipc_weights * (fluc_fit_out$Qn_shift_star - psi)

  # sanity check on EIF
  # NOTE: EIF ~ N(0, V(D(P)(o))), so mean(EIF) ~= 0
  eif_msg <- dplyr::if_else(
    abs(mean(eif)) < tol_eif,
    paste("EIF mean <", tol_eif, "(sufficiently low)."),
    paste(
      "EIF mean =", mean(eif),
      "(higher than expected)."
    )
  )

  # compute the variance based on the EIF
  # NOTE: scale by length of observations to get on same scale as parameter
  # NOTE: var(eif) and mean(eif^2) are nearly equivalent
  # var_eif <- mean(eif^2) / length(Delta)  <= denom since we use full N
  var_eif <- stats::var(eif) / length(Delta)

  # return the variance and the EIF value at each observation
  out <- list(psi = psi, var = var_eif, eif = eif, msg = eif_msg)
  return(out)
}

