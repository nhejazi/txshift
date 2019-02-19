#' Estimate Auxiliary Covariate from Efficient Influence Function
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
  # compute the ratio of the propensity scores for Hn(A,W)
  ratio_g_noshift <- gn$downshift / gn$noshift
  ratio_g_noshift[abs(ratio_g_noshift) == Inf] <- 0
  ratio_g_noshift[is.na(ratio_g_noshift)] <- 0

  # compute the ratio of the propensity scores for Hn(d(A,W),W)
  ratio_g_shift <- gn$noshift / gn$upshift
  ratio_g_shift[abs(ratio_g_shift) == Inf] <- 0
  ratio_g_shift[is.na(ratio_g_shift)] <- 0

  # TODO: consider case where there is not support everywhere
  H_n <- data.table::as.data.table(cbind(ratio_g_noshift, ratio_g_shift))
  data.table::setnames(H_n, c("noshift", "shift"))
  return(H_n)
}
