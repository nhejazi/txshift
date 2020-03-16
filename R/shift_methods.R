#' Simple Additive Modified Treatment Policy
#'
#' @details TODO
#'
#' @param A A \code{numeric} vector of observed treatment values.
#' @param W A \code{numeric} matrix of observed baseline covariate values.
#' @param delta A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the treatment \code{A}.
#'
#' @return TODO
shift_additive <- function(A, W = NULL, delta) {
  shifted_treatment <- A + delta
  return(shifted_treatment)
}
