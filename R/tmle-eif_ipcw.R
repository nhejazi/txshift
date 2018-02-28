#' Compute IPC-Weighted Component of Efficient Influence Function
#'
#' Computes an additional component of the efficient influence function needed
#' for efficient estimation of IPCW-TMLEs. This takes the form:
#' %0=P_n(\Delta-\Pi_n^*(V)))\frac{\E(D^F(P^0_{X,n})\mid\Delta=1,V)}{\Pi^*_n(V)}
#'
#' @param fluc_mod_out ...
#' @param data_in ...
#' @param Hn ...
#' @param Delta ...
#' @param ipc_weights ...
#' @param ipc_weights_norm ...
#' @param eif_tol ...
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
tmle_eif_ipcw <- function(fluc_mod_out,
                          data_in,
                          Hn,
                          Delta,
                          ipc_weights = rep(1, nrow(data_in)),
                          ipc_weights_norm = rep(1, nrow(data_in)),
                          eif_tol = 1e-10) {

  # compute TMLE of the treatment shift parameter
  param_obs_est <- rep(0, length(Delta))
  param_obs_est[Delta == 1] <- ipc_weights_norm * fluc_mod_out$Qn_shift_star
  psi <- sum(param_obs_est)

  # compute the efficient influence function (EIF) / canonical gradient (D*)
  eif <- rep(0, length(Delta))
  eif[Delta == 1] <- ipc_weights * Hn$noshift *
    (data_in$Y - fluc_mod_out$Qn_noshift_star) +
    ipc_weights * (fluc_mod_out$Qn_shift_star - psi)

  # sanity check on EIF
  # NOTE: EIF ~ N(0, V(D(P)(o))), so mean(EIF) ~= 0
  eif_msg <- dplyr::if_else(
    abs(mean(eif)) < eif_tol,
    paste("EIF mean <", eif_tol, "(sufficiently low)."),
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

