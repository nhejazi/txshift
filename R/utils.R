#' Print Method for Counterfactual Mean of Stochastic Shift Intervention
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
#'   g_exp_fit_args = list(
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

###############################################################################

#' Print Method for Marginal Structural Models
#'
#' @details The \code{print} method for objects of class \code{txshift_msm}.
#'
#' @param x An object of class \code{txshift_msm}.
#' @param ... Other options (not currently used).
#'
#' @method print txshift_msm
#'
#' @importFrom tibble as_tibble
#'
#' @return None. Called for the side effect of printing particular slots of
#'  objects of class \code{txshift_msm}.
#'
#' @examples
#' if (require("sl3")) {
#'   set.seed(3287)
#'   n_obs <- 1000
#'   W <- as.numeric(replicate(1, rbinom(n_obs, 1, 0.5)))
#'   A <- as.numeric(rnorm(n_obs, mean = 2 * W, sd = 1))
#'   Y <- rbinom(n_obs, 1, plogis(2 * A - W))
#'   msm <- msm_vimshift(
#'     W = W, A = A, Y = Y, estimator = "tmle",
#'     g_exp_fit_args = list(
#'       fit_type = "sl",
#'       sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
#'     ),
#'     Q_fit_args = list(
#'       fit_type = "glm",
#'       glm_formula = "Y ~ ."
#'     ),
#'     delta_grid = seq(-1, 1, 0.25)
#'   )
#'   print(msm)
#' }
#' @export
print.txshift_msm <- function(x, ...) {
  if (x[["msm_type"]] == "piecewise") {
    print(paste(
      x[["msm_type"]], "MSM with knot point at x =",
      x[["msm_knot"]]
    ))
  } else {
    print(paste(x[["msm_type"]], "MSM"))
  }
  print(tibble::as_tibble(x[["msm_est"]]))
}

###############################################################################

#' Summary for Counterfactual Mean of Stochastic Shift Intervention
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
#'   g_exp_fit_args = list(
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

###############################################################################

is.txshift <- function(x) {
  class(x) == "txshift"
}

is.txshift_msm <- function(x) {
  class(x) == "txshift_msm"
}
