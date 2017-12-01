#' Estimate Inverse Probability of Censoring Weights for IPCW-TMLEs
#'
#' description
#'
#' @param V ...
#' @param Delta ...
#'
#' @importFrom stats glm predict binomial
#'
#' @return A \code{numeric} vector corresponding to the inverse probability of
#'  censoring weights required for computing an IPCW-TMLE. Formally, this is
#'  nothing more than \frac{\Delta}{\Pi_n}, where the term \Pi_n is simply the
#'  predicted probability of belonging to a censoring class as computed using
#'  standard logistic regression.
#
ipcw_est <- funcion(V,
                    Delta) {
  # fit a logistic regression to get class probabilities for IPCW
  ipcw_reg <- stats::glm(Delta ~ V, family = stats::binomial)
  ipcw_probs <- stats::predict(object = ipcw_reg,
                               newdata = as.data.frame(cbind(V, Delta)))

  # compute the inverse weights as Delta/Pi_n and return this vector
  ipc_weights <- (Delta / as.numeric(ipcw_probs))
  ipc_weights_out <- ipc_weights[ipc_weights != 0]
  return(ipc_weights_out)
}

