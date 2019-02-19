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
#' @importFrom stats qnorm
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

  # compute the EIF variance multiplier for the CI
  # NOTE: the variance value is already scaled by length of observations
  sd_eif <- sqrt(object$var)

  # compute the interval around the point estimate
  ci_psi <- norm_bounds * sd_eif + object$psi

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
#'
#' @importFrom stats confint
#'
#' @method summary txshift
#'
#' @export
#
summary.txshift <- function(object,
                            ...,
                            ci_level = 0.95) {

  # compute confidence interval using the pre-defined method
  ci <- stats::confint(object, level = ci_level)

  # only print useful info about the mean of the efficient influence function
  eif_mean <- format(mean(object$eif), scientific = TRUE)

  # create output table from input object and confidence interval results
  out <- c(round(c(ci, object$var), digits = 6), eif_mean, object$n_iter)
  names(out) <- c(
    "lwr_ci", "param_est", "upr_ci", "param_var",
    "eif_mean", "n_iter"
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
  print(x[c("psi", "var", "msg", "n_iter")])
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

#' Scaling by Inducing Boundedness
#'
#' description
#'
#' @param Y A \code{numeric} vector corresponding to the observed values of the
#'  outcome variable of interest.
#' @param pred_vals A \code{numeric} vector corresponding to predicted values of
#'  the outcome of interest (i.e., Qn in the Targeted Learning notation).
#' @param scale_target A \code{numeric} vector specifying the quantity that is
#'  to be re-scaled in the way specified below in \code{scale_type}.
#' @param scale_type An atomic \code{character} vector specifying the type of
#'  scaling to be performed. Use "bound_in_01" to force \code{scale_target}
#'  above to be bounded in the interval (0, 1) with respect to the outcome
#'  \code{Y}. The other option, "observed_vals", re-scales \code{scale_target}
#'  to be on the same scale as the input \code{Y}.
#'
#' @keywords internal
#
bound_scaling <- function(Y,
                          pred_vals = NULL,
                          scale_target = Y,
                          scale_type = c("bound_in_01", "observed_vals")) {
  # check arguments
  scale_type <- match.arg(scale_type)

  # compute minimum and maximum of Y
  y_min <- min(Y)
  y_max <- max(Y)

  if (scale_type == "bound_in_01") {
    out_star <- (scale_target - y_min) / (y_max - y_min)
    return(out_star)
  } else if (scale_type == "observed_vals") {
    out_observed_scale <- (y_max - y_min) * scale_target + y_min
    return(out_observed_scale)
  }
}

