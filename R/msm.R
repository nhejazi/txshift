#' Working marginal structural model for causal effects of an intervention grid
#'
#' @details Computes estimates of the counterfactual mean over a grid of shift
#'  stochastic interventions and fits a working marginal structural model to
#'  summarize the trend through the counterfactual means as a function of the
#'  specified shift intervention. The working marginal structural model may be
#'  linear in the shift parameter or piecewise linear with a single knot point.
#'  Provides support for two weighting schemes, may be used with either of the
#'  one-step or TML estimators, and also allows the construction of marginal or
#'  simultaneous confidence intervals.
#'
#' @param W A \code{matrix}, \code{data.frame}, or similar containing a set of
#'  baseline covariates.
#' @param A A \code{numeric} vector corresponding to a treatment variable. The
#'  parameter of interest is defined as a location shift of this quantity.
#' @param C_cens A \code{numeric} indicator for whether a given observation was
#'  subject to censoring by way of loss to follow-up. The default assumes no
#'  censoring due to loss to follow-up.
#' @param Y A \code{numeric} vector of the observed outcomes.
#' @param C_samp A \code{numeric} indicator for whether a given observation was
#'  subject to censoring by being omitted from the second-stage sample, used to
#'  compute an inverse probability of censoring weighted estimator in such
#'  cases. The default assumes no censoring due to two-phase sampling.
#' @param V The covariates that are used in determining the sampling procedure
#'  that gives rise to censoring. The default is \code{NULL} and corresponds to
#'  scenarios in which there is no censoring (in which case all values in the
#'  preceding argument \code{C} must be uniquely 1. To specify this, pass in a
#'  NAMED \code{list} identifying variables amongst W, A, Y that are thought to
#'  have played a role in defining the sampling/censoring mechanism (C).
#' @param delta_grid A \code{numeric} vector giving the individual values of
#'  the shift parameter used in computing each of the estimates.
#' @param msm_form A \code{list} specifying the type of working MSM to fit to
#'  summarize the counterfactual means. The \code{list} has two components:
#'  (1) \code{"type"}, which may be either "linear" or "piecewise", and (2)
#'  \code{"knot"}, which, if specified, must be a value in \code{delta_grid}.
#'  See examples for its use.
#' @param estimator The type of estimator to be fit, either \code{"tmle"} for
#'  targeted maximum likelihood estimation or \code{"onestep"} for a one-step
#'  augmented inverse probability weighted (AIPW) estimator.
#' @param weighting Whether to weight each parameter estimate by the inverse of
#'  its variance (in order to improve stability of the resultant MSM fit) or to
#'  simply weight all parameter estimates equally. The default is the option
#'  \code{"identity"}, weighting all estimates identically.
#' @param ci_level A \code{numeric} indicating the desired coverage level of
#'  the confidence interval to be computed.
#' @param ci_type Whether to construct a simultaneous confidence band covering
#'  all parameter estimates at once or marginal confidence intervals covering
#'  each parameter estimate separately. The default is to construct marginal
#'  confidence intervals for each parameter estimate rather than a simultaneous
#'  confidence band.
#' @param ... Additional arguments to be passed to \code{\link{txshift}}.
#'
#' @importFrom assertthat assert_that
#' @importFrom data.table as.data.table copy setnames
#' @importFrom stats cov confint qnorm pnorm as.formula model.matrix lm
#' @importFrom lspline lspline
#' @importFrom mvtnorm qmvnorm
#'
#' @return A \code{list} containing estimates of the individual counterfactual
#'  means over a grid in the shift parameters (\code{delta_grid}), alongside
#'  the estimate of a marginal structural model that summarizes a trend through
#'  these counterfactual means.
#'
#' @examples
#' if (require("sl3")) {
#'   n_obs <- 100
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
#'
#'   # fit a linear spline with knot at 0
#'   n_obs <- 100
#'   W <- as.numeric(replicate(1, rbinom(n_obs, 1, 0.5)))
#'   A <- as.numeric(rnorm(n_obs, mean = 2 * W, sd = 1))
#'   Y <- rbinom(n_obs, 1, plogis(0.1 * A * (A >= 0) - 3 * A * (A < 0) - W))
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
#'     delta_grid = seq(-1, 1, 0.25),
#'     msm_form = list(type = "piecewise", knot = 0)
#'   )
#' }
#' @export
msm_vimshift <- function(W,
                         A,
                         C_cens = rep(1, length(Y)),
                         Y,
                         C_samp = rep(1, length(Y)),
                         V = NULL,
                         delta_grid = seq(-0.5, 0.5, 0.5),
                         msm_form = list(type = "linear", knot = NA),
                         estimator = c("tmle", "onestep"),
                         weighting = c("identity", "variance"),
                         ci_level = 0.95,
                         ci_type = c("marginal", "simultaneous"),
                         ...) {
  # set default values of arguments
  estimator <- match.arg(estimator)
  weighting <- match.arg(weighting)
  ci_type <- match.arg(ci_type)

  # sanity checks for MSM
  assertthat::assert_that(msm_form[["type"]] %in% c("linear", "piecewise"))
  if (msm_form[["type"]] == "piecewise") {
    assertthat::assert_that(!is.na(msm_form[["knot"]]),
      msg = "knot cannot be NA for piecewise MSMs"
    )
  }

  # make sure more than one parameter is to be estimated for trend test
  assertthat::assert_that(length(delta_grid) > 1)

  # fit TML or one-step estimator for each value of shift in the grid
  est_over_grid <-
    lapply(delta_grid, function(shift) {
      est <- txshift(
        W = W, A = A, Y = Y, C_samp = C_samp, V = V,
        delta = shift, estimator = estimator,
        ...
      )
    })

  # matrix of EIF(O_i) values and estimates across each parameter estimated
  eif_mat <- do.call(cbind, lapply(est_over_grid, `[[`, "eif"))
  psi_vec <- do.call(c, lapply(est_over_grid, `[[`, "psi"))

  # multiplier for CI construction: simultaneous confidence interval
  if (ci_type == "simultaneous" && (ncol(eif_mat) > 1)) {
    # compute correlation based on covariance of EIF
    var_eif <- stats::cov(eif_mat)
    rho_eif <- var_eif / sqrt(tcrossprod(diag(var_eif)))
    mvt_eif <- mvtnorm::qmvnorm(ci_level, tail = "both", corr = rho_eif)
    # NOTE: c(-1, 1) instead of c(1, -1): mvtnorm call differs from qnorm
    ci_mult <- c(-1, 1) * mvt_eif$quantile
  } else {
    ci_mult <- c(1, -1) * stats::qnorm((1 - ci_level) / 2)
  }
  # create confidence intervals, overriding default multiplier
  wald_cis <- do.call(rbind, lapply(est_over_grid, function(est) {
    stats::confint(est, ci_mult = ci_mult)
  }))

  # set weights to be the inverse of the variance of each TML estimate
  if (weighting == "variance") {
    weights <- as.numeric(1 / diag(stats::cov(eif_mat)))
  } else {
    weights <- rep(1, length(psi_vec))
  }

  # set right-hand side of MSM formula
  if (msm_form[["type"]] == "piecewise") {
    msm_rhs <- paste0("delta + I(pmax(delta - ", msm_form[["knot"]], ", 0))")
  } else {
    msm_rhs <- "delta"
  }

  # create design matrix for MSM
  x_mat <- stats::model.matrix(
    stats::as.formula(paste0("psi_vec ~ ", msm_rhs)),
    data = data.frame(psi_vec = psi_vec, delta = delta_grid)
  )

  # compute the MSM parameters
  omega <- diag(weights)
  s_mat <- solve(t(x_mat) %*% omega %*% x_mat) %*% t(x_mat) %*% omega
  msm_param <- as.vector(s_mat %*% psi_vec)

  # compute inference for MSM based on individual EIF(O_i) for each parameter
  msm_eif <- tcrossprod(eif_mat, s_mat)
  msm_var <- diag(stats::cov(msm_eif))
  msm_se <- sqrt(msm_var / nrow(msm_eif))

  # build confidence intervals and hypothesis tests for EIF(msm)
  ci_msm_param <- msm_se %*% t(ci_mult) + msm_param
  pval_msm_param <- 2 * stats::pnorm(-abs(msm_param / msm_se))

  # summarize output of individual shift-specific estimates
  vimshift_out <- data.table::as.data.table(
    list(
      delta = delta_grid,
      ci_lwr = wald_cis[, 1],
      psi = psi_vec,
      ci_upr = wald_cis[, 3]
    )
  )

  # create summary table for MSM estimates
  msm_out <- data.table::as.data.table(list(
    param = names(msm_se),
    ci_lwr = ci_msm_param[, 1],
    param_est = msm_param,
    ci_upr = ci_msm_param[, 2],
    param_se = msm_se,
    p_value = pval_msm_param
  ))

  # create and rename MSM data for downstream ggplot2 compatibility
  msm_data <- data.table::copy(vimshift_out)
  data.table::setnames(msm_data, c("psi", "delta"), c("y", "x"))

  # compute linear working MSM or single-knot spline model
  if (msm_form[["type"]] == "piecewise" && !is.na(msm_form[["knot"]])) {
    msm_fit <- stats::lm(y ~ lspline::lspline(x, msm_form[["knot"]],
      marginal = TRUE
    ),
    weights = weights, data = msm_data
    )
  } else if (msm_form[["type"]] == "linear") {
    msm_fit <- stats::lm(y ~ x, weights = weights, data = msm_data)
  }

  # complete output for MSM
  out <- list(
    param_est = vimshift_out,
    msm_est = msm_out,
    msm_type = msm_form[["type"]],
    msm_data = msm_data,
    msm_fit = msm_fit,
    estimator = estimator,
    delta_grid = delta_grid,
    ci_type = ci_type,
    ci_level = ci_level
  )
  class(out) <- "txshift_msm"
  return(out)
}
