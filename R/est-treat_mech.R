#' Estimate the Treatment Mechanism
#'
#' Compute the propensity score (treatment mechanism) for the observed data,
#' including the shift. This returns the propensity score for the observed data
#' (at A_i) and the shift of interest (at A_i - delta).
#'
#' @param A A \code{numeric} vector of observed treatment values.
#' @param W A \code{numeric} matrix of observed baseline covariate values.
#'
#' @importFrom condensier fit_density predict_probability
#'
#' @keywords internal
#'
#' @export
#'
#' @author Nima Hejazi
#' @author David Benkeser
#
est_g <- function(A, W) {
   data_O <- as.data.frame(cbind(A, W))
   dens_fit <- condensier::fit_density(X = paste0("W", seq_len(W)), Y = "A",
                                       input_data = data_O, nbins = 20,
                                       bin_method = "equal.mass",
                                       bin_estimator = condensier::speedglmR6$new()
                                      )
   #...
}

