#' Estimate the Auxiliary ("Clever") Covariate
#'
#' description
#'
#' @param gn ...
#' @param a ...
#' @param w ...
#'
#' @keywords internal
#
est_Hn <- function(gn, a = NULL, w = NULL) {
  # compute upper and lower limits for treatment
  # ...
  # ...

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
  H_n <- as.data.frame(cbind(ratio_g_noshift, ratio_g_shift))
  colnames(H_n) <- c("noshift", "shift")

  # output
  return(H_n)
}
