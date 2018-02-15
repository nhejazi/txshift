utils::globalVariables(c("key"))

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

  # find duplicated rows of data_in (removing "Y" column)
  rm_col <- "Y"
  dups <- duplicated(data_in[ , !rm_col, by = key(data_in), with = FALSE])
  duplicated_idx <- which(dups)

  # if no duplicates duplicated_idx == integer(0)
  if (length(duplicated_idx) == 0) {
    psi <- sum(fluc_fit_out$Qn_shift_star * eps_updated$qn[ipc_weights > 0])
  } else {
    psi <- sum(fluc_fit_out$Qn_shift_star[-duplicated_idx] *
               eps_updated$qn[-duplicated_idx])
  }

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

