#' Compute IPCW Part of the Efficient Influence Function for IPCW-TMLE
#'
#' Computes an additional component of the efficient influence function needed
#' for efficient estimation of IPCW-TMLEs. This takes the form:
#' %0=P_n(\Delta-\Pi_n^*(V)))\frac{\E(D^F(P^0_{X,n})\mid\Delta=1,V)}{\Pi^*_n(V)}
#'
#' @param fluc_fit_out ...
#' @param eps_updated ...
#' @param data_in ...
#' @param Hn ...
#' @param Y A \code{numeric} vector of observed outcomes.
#' @param Delta ...
#' @param ipc_weights ...
#' @param tol_eif ...
#'
#' @importFrom stats var
#' @importFrom dplyr if_else add_count select mutate
#' @importFrom data.table is.data.table as.data.table set
#'
#' @keywords internal
#'
#' @export
#'
#' @author Nima Hejazi
#' @author David Benkeser
#
tmle_eif_ipcw <- function(fluc_fit_out,
                          eps_updated,
                          data_in,
                          Hn,
                          Y,
                          Delta,
                          ipc_weights = rep(1, length(Y)),
                          tol_eif = 1e-7) {
  # create data set to facilitate computing the TMLE for Psi
  if (!data.table::is.data.table(data_in)) {
    data_in <- data.table::as.data.table(data_in)
  }
  data.table::set(data_in, j = "Qn", value = fluc_fit_out$Qn_shift_star)
  data.table::set(data_in, j = "qn",
                  value = eps_updated$qn[eps_updated$qn != 0])

  # find the contributions to Psi from the unique observations
  dat_with_psi <- dat %>%
    dplyr::add_count(W, round(A, 2)) %>%
    dplyr::select(Qn, qn, n) %>%
    dplyr::mutate(
      # let's divide by n here since, on average, this'll catch the contribution
      # as though the observation were unique (???).
      # note: doesn't affect unique contributions because division by 1
      psi_est = (Qn * qn) / n
    )
  # ...
  psi_obs_est[Delta == 1] <- dat_with_psi$psi_est
  psi <- sum(dat_with_psi$psi_est)

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

