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

################################################################################

#' Summary for Counterfactual Mean Under Stochastic Intervention
#'
#' @details Print a convenient summary for objects computed using
#'  \code{\link{txshift}}.
#'
#' @param object An object of class \code{txshift}, as produced by invoking
#'  the function \code{\link{txshift}}, for which a confidence interval is to
#'  be computed.
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
#' @return None. Called for the side effect of printing a summary of particular
#'  slots of objects of class \code{txshift}.
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
#' summary(txout)
#' @export
summary.txshift <- function(object,
                            ...,
                            ci_level = 0.95,
                            digits = 4) {
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
#' @details The \code{print} method for objects of class \code{txshift}.
#'
#' @param x An object of class \code{txshift}.
#' @param ... Other options (not currently used).
#'
#' @method print txshift
#'
#' @return None. Called for the side effect of printing particular slots of
#'  objects of class \code{txshift}.
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
#' print(txout)
#' @export
print.txshift <- function(x, ...) {
  print(x[c("psi", "var", "estimator", "n_iter")])
}

################################################################################

#' Bound Precision
#'
#' @details Bound values in the unit interval to machine precision in order to
#'  avoid numerical instability issues in downstream computation.
#'
#' @param vals \code{numeric} vector of values in the interval [0, 1] to be
#'  bounded within arbitrary machine precision. The most common use of this
#'  functionality is to avoid indeterminate or non-finite values after the
#'  application \code{stats::qlogis}.
#'
#' @importFrom assertthat assert_that
#'
#' @return A \code{numeric} vector of the same length as \code{vals}, where
#'  the returned values are bounded to machine precision. This is intended to
#'  avoid numerical instability issues.
bound_precision <- function(vals) {
  assertthat::assert_that(!(max(vals) > 1 | min(vals) < 0))
  vals[vals == 0] <- .Machine$double.neg.eps
  vals[vals == 1] <- 1 - .Machine$double.neg.eps
  return(vals)
}

################################################################################

#' Bound Generalized Propensity Score
#'
#' @details Bound estimated values of the generalized propensity score (a
#'  conditional density) to avoid numerical instability issues arising from
#'  practical violations of the assumption of positivity.
#'
#' @param vals \code{numeric} vector of propensity score estimate values. Note
#'  that, for this parameter, the propensity score is (conditional) density and
#'  so it ought not be bounded from above.
#'
#' @return A \code{numeric} vector of the same length as \code{vals}, where the
#'  returned values are bounded such that the minimum is no lower than 1/n, for
#'  the sample size n.
bound_propensity <- function(vals) {
  # bound likelihood component g(a|w) away from 0 only
  propensity_bound <- 1 / length(vals)
  vals[vals < propensity_bound] <- propensity_bound
  return(vals)
}

################################################################################

#' Transform values by scaling to the unit interval
#'
#' @details A transformation that scales an arbitrary set of input values to
#'  the unit interval. See \code{\link{scale_to_original}} for a corresponding
#'  backtransformation.
#'
#' @param vals A \code{numeric} vector corresponding to the observed values of
#'  the variable of interest, to be re-scaled to the unit interval [0,1].
#'
#' @return A \code{numeric} vector of the same length as \code{vals}, where the
#'  values are re-scaled to lie in unit interval [0, 1].
scale_to_unit <- function(vals) {
  # compute re-scaled value in interval [0,1]
  scaled_vals <- (vals - min(vals)) / (max(vals) - min(vals))
  return(scaled_vals)
}

################################################################################

#' Transform values from the unit interval back to their original scale
#'
#' @details A back-transformation that returns values computed in the unit
#'  interval to their original scale. This is used in re-scaling updated TML
#'  estimates back to their natural scale. Undoes \code{\link{scale_to_unit}}.
#'
#' @param scaled_vals A \code{numeric} vector corresponding to re-scaled values
#'  in the unit interval, to be re-scaled to the original interval.
#' @param max_orig A \code{numeric} scalar value giving the maximum of the
#'  values on the original scale.
#' @param min_orig A \code{numeric} scalar value giving the minimum of the
#'  values on the original scale.
#'
#' @return A \code{numeric} vector of the same length as \code{scaled_vals},
#'  where the values are re-scaled to lie in their original/natural interval.
scale_to_original <- function(scaled_vals, max_orig, min_orig) {
  scaled_orig <- scaled_vals * (max_orig - min_orig) + min_orig
  return(scaled_orig)
}
