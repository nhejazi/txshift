#' Confidence Intervals for Counterfactual Mean Under Stochastic Intervention
#'
#' @details Compute confidence intervals for estimates produced by
#'  \code{\link{txshift}}.
#'
#' @param object An object of class \code{txshift}, as produced by invoking
#'  the function \code{\link{txshift}}, for which a confidence interval is to
#'  be computed.
#' @param parm A \code{numeric} vector indicating indices of \code{object$est}
#'  for which to return confidence intervals.
#' @param level A \code{numeric} indicating the level of the confidence
#'  interval to be computed.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom stats qnorm plogis qlogis
#'
#' @method confint txshift
#'
#' @return A named \code{numeric} vector containing the parameter estimate from
#'  a \code{txshift} object, alongside lower and upper Wald-style confidence
#'  intervals at a specified coverage level.
#'
#' @examples
#' set.seed(429153)
#' n_obs <- 100
#' W <- replicate(2, rbinom(n_obs, 1, 0.5))
#' A <- rnorm(n_obs, mean = 2 * W, sd = 1)
#' Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
#' txout <- txshift(
#'   W = W, A = A, Y = Y, delta = 0.5,
#'   estimator = "tmle",
#'   g_fit_args = list(
#'     fit_type = "hal", n_bins = 5,
#'     grid_type = "equal_mass",
#'     lambda_seq = exp(-1:-9)
#'   ),
#'   Q_fit_args = list(
#'     fit_type = "glm",
#'     glm_formula = "Y ~ ."
#'   )
#' )
#' confint(txout)
#' @export
confint.txshift <- function(object,
                            parm = seq_len(object$psi),
                            level = 0.95,
                            ...) {
  # first, let's get Z_(1 - alpha)
  norm_bounds <- c(-1, 1) * abs(stats::qnorm(p = (1 - level) / 2))

  if (length(unique(object$outcome)) > 2) { # assume continuous outcome
    # compute the EIF variance multiplier for the CI
    # NOTE: the variance value is already scaled by length of observations
    sd_eif <- sqrt(object$var)

    # compute the interval around the point estimate
    ci_psi <- norm_bounds * sd_eif + object$psi
  } else if (length(unique(object$outcome)) == 2) { # binary outcome
    # for binary outcome case, compute on the logit scale and back-transform
    psi_ratio <- stats::qlogis(object$psi)
    grad_ratio_delta <- (1 / object$psi) + (1 / (1 - object$psi))
    se_eif_logit <- sqrt(grad_ratio_delta^2 * object$var)
    ci_psi <- stats::plogis(norm_bounds * se_eif_logit + psi_ratio)
  } else {
    stop("The outcome has fewer than 2 levels: this case is not supported.")
  }

  # set up output CI object
  ci_out <- c(ci_psi[1], object$psi, ci_psi[2])
  names(ci_out) <- c("lwr_ci", "est", "upr_ci")
  return(ci_out)
}
