#' Simple Additive Modified Treatment Policy
#'
#' @details A simple modified treatment policy that modifes the observed value
#'  of the exposure by shifting it by a value \code{delta}. Note that this
#'  shifting function assumes support of A|W across all strata of W.
#'
#' @param A A \code{numeric} vector of observed treatment values.
#' @param W A \code{numeric} matrix of observed baseline covariate values.
#' @param delta A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the treatment \code{A}.
#'
#' @return A \code{numeric} vector containing the shifted exposure values.
shift_additive <- function(A, minA, maxA, W = NULL, delta) {
  if(delta < 0){
    shifted_treatment <- ifelse((minA<= A + delta),  A + delta, A)
  } else {
    shifted_treatment <- ifelse((maxA>= A + delta),  A + delta, A)
  }
 
  return(shifted_treatment)
}
