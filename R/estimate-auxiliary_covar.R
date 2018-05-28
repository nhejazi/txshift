#' Estimate the Auxiliary ("Clever") Covariate
#'
#' @param gn An estimate of the treatment probability (propensity score), using
#'  the output provided by internal function \code{estimate-propensity_score}.
#' @param a A \code{numeric} vector of the observed values of the treatment.
#' @param w A \code{numeric}, \code{matrix}, \code{data.frame} or similar object
#'  that contains the observed values of the baseline covariates.
#'
#' @importFrom data.table as.data.table setnames
#'
#' @export
#
est_Hn <- function(gn, a = NULL, w = NULL) {
  # TODO: compute upper and lower limits for treatment

  # compute the ratio of the propensity scores for Hn(A,W)
  ratio_g_noshift <- gn$downshift / gn$noshift

  # compute the ratio of the propensity scores for Hn(d(A,W),W)
  ratio_g_shift <- gn$noshift / gn$upshift

  # modify the ratio of the propensity scores
  # based on the indicators for shifting
  # ind_a <- ...
  # ind_a_delta <- ...
  # H_n <- ind_a * ratio_g + ind_a_delta

  # TODO: consider case where there is not support everywhere
  # that is, when the indicators kick in -- ignored for now...
  H_n <- data.table::as.data.table(cbind(ratio_g_noshift, ratio_g_shift))
  data.table::setnames(H_n, c("noshift", "shift"))

  # output
  return(H_n)
}
