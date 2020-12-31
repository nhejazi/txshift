#' Print Method for Counterfactual Mean of Stochastic Shift Intervention
#'
#' @details The \code{print} method for objects of class \code{txshift}.
#'
#' @param x An object of class \code{txshift}.
#' @param ... Other options (not currently used).
#' @param ci_level A \code{numeric} indicating the level of the confidence
#'  interval to be computed.
#'
#' @method print txshift
#'
#' @importFrom stats confint
#' @importFrom scales percent
#' @importFrom cli cli_text col_red
#'
#' @return None. Called for the side effect of printing an informative summary
#'  of slots of objects of class \code{txshift}.
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
print.txshift <- function(x, ..., ci_level = 0.95) {
  # compute confidence interval
  ci <- stats::confint(x, level = ci_level)

  # construct and print output
  cli::cli_text("{.strong Counterfactual Mean for Shifted Treatment}")
  cli::cli_text(cat("    "), "{.strong Intervention}: Treatment + ",
                cli::col_red("{x$delta}"))
  cat("\n")
  cli::cli_text("{.strong txshift Estimator}: {x$estimator}")
  cli::cli_text(cat("    "), "{.strong Estimate}: {round(x$psi, 4)}")
  cli::cli_text(cat("    "), "{.strong Std. Error}: {round(sqrt(x$var), 4)}")
  cli::cli_text(cat("    "), "{.strong {scales::percent(ci_level)} CI}:
                [{round(ci[1], 4)}, {round(ci[3], 4)}]")
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
#' @importFrom cli cli_text col_red col_blue
#'
#' @return None. Called for the side effect of printing an informative summary
#'  of slots of objects of class \code{txshift_msm}.
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
  # construct and print output
  cli::cli_text("{.strong {x$msm_type} MSM over Grid of Shifts in Treatment}")
  cli::cli_text(cat("    "), "{.strong Intervention Grid}: Treatment + (",
                cli::col_red("{x$delta_grid}"), ")")
  if (x[["msm_type"]] == "piecewise") {
    cli::cli_text(cat("    "), "{.strong Knot point}: ",
                  cli::col_blue("{x$msm_knot}"))
  }
  cat("\n")
  cli::cli_text("{.strong txshift MSM Estimator}: {x$estimator}")
  cli::cli_text(cat("    "), "{.strong Estimated Slope}:
                {round(x$msm_est$param_est[2], 4)}")
  cli::cli_text(cat("    "), "{.strong Std. Error}:
                {round(x$msm_est$param_se[2], 4)}")
  cli::cli_text(cat("    "), "{.strong {scales::percent(x$ci_level)} CI}:
                [{round(x$msm_est$ci_lwr[2], 4)},
                 {round(x$msm_est$ci_upr[2], 4)}]")
  cli::cli_text(cat("    "), "{.strong p-value (vs. no trend)}:
                {round(x$msm_est$p_value[2], 4)}")
}

###############################################################################

is.txshift <- function(x) {
  class(x) == "txshift"
}

is.txshift_msm <- function(x) {
  class(x) == "txshift_msm"
}
