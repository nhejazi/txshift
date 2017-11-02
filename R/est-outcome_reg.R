#' Estimate the Outcome Regression
#'
#' Compute the outcome regression for the observed data, including with the
#' shift imposed by the intervention. This returns the propensity score for the
#' observed data (at A_i) and the shift of interest (at A_i - delta).
#'
#' @param Y A \code{numeric} vector of observed outcomes.
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
est_Q <- function(Y, A, W, fit_type = c("glm", "sl"), outcome_type) {
   # generate a single data structure from the inputs
   data_O <- as.data.frame(cbind(Y, A, W))

   # compute the outcome regression over the observed data
   if (fit_type == "glm") {
      mod_Q <- stats::glm(Y ~ A + W, data = data_O)
   } else if (fit_type == "sl") {
      mod_Q <- SuperLearner::SuperLearner(Y ~ A + W, data = data_O)
   }

   # get Qn(a,w)
   out_QnAW <- predict(mod_Q, newdata = data_O)

   # shift the treatment values
   data_O_shift <- data_O
   data_O_shift$A <- data_O$A + delta

   # compute the outcome regression under the treatment shift
   out_QndAW <- predict(mod_Q, newdata = data_O_shift)

   # creat output matrix
   out <- as.data.frame(cbind(out_QnAW, outQndAW))
   return(out)
}

