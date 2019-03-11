#' Simple Modified Treatment Policy (Shifted Treatment)
#'
#' @param A A \code{numeric} vector of observed treatment values.
#' @param W A \code{numeric} matrix of observed baseline covariate values.
#' @param delta A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the treatment \code{A}.
#'
#' @keywords internal
#
tx_shift <- function(A,
                     W = NULL,
                     delta) {
  # could support more than just additive shifts?
  shifted_A <- A + delta
  return(shifted_A)
}

################################################################################

#' Confidence Intervals for Shifted Treatment Parameters
#'
#' Compute confidence intervals for estimates produced by \code{tmle_txshift}
#'
#' @param object An object of class \code{txshift}, as produced by invoking
#'  the function \code{tmle_txshift}, for which a confidence interval is to be
#'  computed.
#' @param parm A \code{numeric} vector indicating indices of \code{object$est}
#'  for which to return confidence intervals.
#' @param level A \code{numeric} indicating the level of the confidence interval
#'  to be computed.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom stats qnorm plogis qlogis
#'
#' @method confint txshift
#'
#' @export
#
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

################################################################################

#' Summary for Shifted Treatment Parameter Objects
#'
#' Print a convenient summary for objects computed using \code{tmle_txshift}.
#'
#' @param object An object of class \code{txshift}, as produced by invoking
#'  the function \code{tmle_txshift}, for which a confidence interval is to be
#'  computed.
#' @param ... Other arguments. Not currently used.
#' @param ci_level A \code{numeric} indicating the level of the confidence
#'  interval to be computed.
#' @param digits A \code{numeric} scalar giving the number of digits to be
#'  displayed or to round results to.
#'
#' @importFrom stats confint
#'
#' @method summary txshift
#'
#' @export
#
summary.txshift <- function(object,
                            ...,
                            ci_level = 0.95,
                            digits = 5) {

  # compute confidence interval using the pre-defined method
  ci <- stats::confint(object, level = ci_level)

  # only print useful info about the mean of the efficient influence function
  eif_mean <- formatC(mean(object$eif), digits = digits, format = "e")

  # create output table from input object and confidence interval results
  out <- c(
    round(c(ci, object$var), digits = digits), eif_mean,
    object$estimator, object$n_iter
  )
  names(out) <- c(
    "lwr_ci", "param_est", "upr_ci", "param_var",
    "eif_mean", "estimator", "n_iter"
  )
  print(noquote(out))
}

################################################################################

#' Print Method for txshift Objects
#'
#' The \code{print} method for objects of class \code{txshift}.
#'
#' @param x An object of class \code{txshift}.
#' @param ... Other options (not currently used).
#'
#' @export
#'
#' @method print txshift
#'
#
print.txshift <- function(x, ...) {
  print(x[c("psi", "var", "estimator", "msg", "n_iter")])
}

################################################################################

#' Bound Precision
#'
#' description
#'
#' @param vals \code{numeric} vector of values in the interval [0, 1] to be
#'  bounded within arbitrary machine precision. The most common use of this
#'  functionality is to avoid indeterminate or non-finite values after the
#'  application \code{stats::qlogis}.
#'
#' @importFrom assertthat assert_that
#'
#' @keywords internal
#
bound_precision <- function(vals) {
  assertthat::assert_that(!(max(vals) > 1 | min(vals) < 0))
  vals[vals == 0] <- .Machine$double.neg.eps
  vals[vals == 1] <- 1 - .Machine$double.neg.eps
  return(vals)
}

################################################################################

#' Scale values to the unit interval [0, 1]
#'
#' @param vals A \code{numeric} vector corresponding to the observed values of
#'  the variable of interest, to be re-scaled to the unit interval [0, 1].
#'
#' @keywords internal
#
scale_to_unit <- function(vals) {
  # compute re-scaled value in interval [0,1]
  scaled_vals <- (vals - min(vals)) / (max(vals) - min(vals))
  return(scaled_vals)
}

################################################################################

#' Scale transformed values in the unit interval to the original scale
#'
#' @param scaled_vals A \code{numeric} vector corresponding to re-scaled values
#'  in the unit interval, to be re-scaled to the original interval.
#' @param max_orig A \code{numeric} scalar value giving the maximum of the
#'  values on the original scale.
#' @param min_orig A \code{numeric} scalar value giving the minimum of the
#'  values on the original scale.
#'
#' @keywords internal
#
scale_to_original <- function(scaled_vals, max_orig, min_orig) {
  scaled_orig <- scaled_vals * (max_orig - min_orig) + min_orig
  return(scaled_orig)
}
