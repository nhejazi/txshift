#' Compute the Shift Parameter Estimate and the Efficient Influence Function
#'
#' Computes the value of the treatment shift parameter as well as statistical
#' inference for the parameter based on the efficient influence function of the
#' target parameter, which takes the following form:
#' %D(P)(o) = H(a,w)(y - \bar{Q}(a,w)) + \bar{Q}(d(a,w)) - \psi(P)
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
tmle_eif <- function(fluc_fit_out,
                     Hn,
                     Y,
                     ipc_weights = rep(1, length(Y)),
                     tol_eif = 1e-7) {
  # compute TMLE
  psi <- mean(ipc_weights * fluc_fit_out$Qn_shift_star)

  # compute the efficient influence function (EIF) / canonical gradient (D*)
  eif <- ipc_weights * Hn$noshift * (Y - fluc_fit_out$Qn_noshift_star) +
    fluc_fit_out$Qn_shift_star - psi

  # sanity check on EIF
  # NOTE: EIF ~ N(0, V(D(P)(o))), so mean(EIF) ~= 0
  eif_msg <- dplyr::if_else(mean(eif) < tol_eif,
                            paste("EIF mean <", tol_eif, "(sufficiently low)."),
                            paste("EIF mean =", mean(eif),
                                  "(higher than expected).")
                           )

  # compute the variance based on the EIF
  # NOTE: scale by length of observations to get on same scale as parameter
  # NOTE: var(eif) and mean(eif^2) are nearly equivalent
  # var_eif <- mean(eif^2) / length(Y)
  var_eif <- stats::var(eif) / length(Y)

  # return the variance and the EIF value at each observation
  out <- list(psi = psi, var = var_eif, eif = eif, msg = eif_msg)
  return(out)
}

