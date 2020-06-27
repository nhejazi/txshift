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
#' @param Y A \code{numeric} vector of the observed outcomes.
#' @param A A \code{numeric} vector corresponding to a treatment variable. The
#'  parameter of interest is defined as a location shift of this quantity.
#' @param W A \code{matrix}, \code{data.frame}, or similar corresponding to a
#'  set of baseline covariates.
#' @param C A \code{numeric} indicator for whether a given observation was
#'  subject to censoring in the two-phase sample. This is used to compute an
#'  IPCW-TMLE in such cases. The default assumes no censoring.
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
#' @param ... Additional arguments to be passed to \code{txshift}.
#'
#' @importFrom assertthat assert_that
#' @importFrom data.table as.data.table
#' @importFrom stats cov confint qnorm pnorm as.formula model.matrix
#' @importFrom mvtnorm qmvnorm
#'
#' @return A \code{list} containing estimates of the individual counterfactual
#'  means over a grid in the shift parameters (\code{delta_grid}), alongside
#'  the estimate of a marginal structural model that summarizes a trend through
#'  these counterfactual means.
#'
#' @examples
#' library(sl3)
#' n_obs <- 100
#' W <- as.numeric(replicate(1, rbinom(n_obs, 1, 0.5)))
#' A <- as.numeric(rnorm(n_obs, mean = 2 * W, sd = 1))
#' Y <- rbinom(n_obs, 1, plogis(2 * A - W))
#' msm <- msm_vimshift(
#'   W = W, A = A, Y = Y, estimator = "tmle",
#'   g_fit_args = list(
#'     fit_type = "sl",
#'     sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
#'   ),
#'   Q_fit_args = list(
#'     fit_type = "glm",
#'     glm_formula = "Y ~ ."
#'   ),
#'   delta_grid = seq(-1, 1, 0.25)
#' )
#'
#' # fit a linear spline with knot at 0
#' n_obs <- 100
#' W <- as.numeric(replicate(1, rbinom(n_obs, 1, 0.5)))
#' A <- as.numeric(rnorm(n_obs, mean = 2 * W, sd = 1))
#' Y <- rbinom(n_obs, 1, plogis(0.1 * A * (A >= 0) - 3 * A * (A < 0) - W))
#' msm <- msm_vimshift(
#'   W = W, A = A, Y = Y, estimator = "tmle",
#'   g_fit_args = list(
#'     fit_type = "sl",
#'     sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
#'   ),
#'   Q_fit_args = list(
#'     fit_type = "glm",
#'     glm_formula = "Y ~ ."
#'   ),
#'   delta_grid = seq(-1, 1, 0.25),
#'   msm_form = list(type = "piecewise", knot = 0)
#' )
#' @export
msm_vimshift <- function(Y,
                         A,
                         W,
                         C = rep(1, length(Y)),
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

  # multiplier for CI construction
  ci_mult <- c(1, -1) * stats::qnorm((1 - ci_level) / 2)

  # fit TML or one-step estimator for each value of shift in the grid
  est_over_grid <-
    lapply(delta_grid, function(shift) {
      est <- txshift(
        W = W, A = A, Y = Y, C = C, V = V,
        delta = shift, estimator = estimator,
        ...
      )
      est_with_ci <- stats::confint(est)
      return(list(
        est_with_ci = est_with_ci,
        eif_from_est = est[["eif"]]
      ))
    })

  # matrix of EIF(O_i) values and estimates across each parameter estimated
  eif_mat <- do.call(cbind, lapply(est_over_grid, `[[`, "eif_from_est"))
  psi_with_ci <- do.call(rbind, lapply(est_over_grid, `[[`, "est_with_ci"))
  psi_vec <- psi_with_ci[, 2]

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

  # simultaneous confidence interval
  if (ci_type == "simultaneous" && (ncol(eif_mat) > 1)) {
    # compute correlation based on covariance of EIF
    var_eif <- stats::cov(eif_mat)
    rho_eif <- var_eif / sqrt(tcrossprod(diag(var_eif)))
    mvt_eif <- mvtnorm::qmvnorm(ci_level, tail = "both", corr = rho_eif)
    # for simultaneous interval, update the quantiles for the CI
    # NOTE: c(-1, 1) instead of c(1, -1): mvtnorm call differs from qnorm
    ci_mult <- c(-1, 1) * mvt_eif$quantile
  }

  # build confidence intervals and hypothesis tests for EIF(msm)
  ci_msm_param <- msm_se %*% t(ci_mult) + msm_param
  pval_msm_param <- 2 * stats::pnorm(-abs(msm_param / msm_se))

  # summarize output of individual shift-specific estimates
  vimshift_out <- data.table::as.data.table(
    list(
      delta = delta_grid,
      ci_lwr = psi_with_ci[, 1],
      psi = psi_vec,
      ci_upr = psi_with_ci[, 3]
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

  # complete output for MSM
  out <- list(
    param_est = vimshift_out,
    msm_est = msm_out,
    msm_type = msm_form[["type"]],
    msm_knot = msm_form[["knot"]]
  )
  class(out) <- "txshift_msm"
  return(out)
}

###############################################################################

#' Plot working MSM for causal effects of an intervention grid
#'
#' @details Creates a visualization of the intervention-specific counterfactual
#'  means as well as the working marginal structural model summarizing the
#'  trend across posited values of the intervention.
#'
#' @param x Object of class \code{txshift_msm} as produced by a call to
#'  \code{\link{msm_vimshift}}.
#' @param ... Additional arguments passed to \code{plot} as necessary.
#'
#' @importFrom ggplot2 ggplot geom_point geom_errorbar geom_abline geom_smooth
#'  aes_string labs theme_bw
#' @importFrom data.table copy setnames
#' @importFrom stats formula lm
#' @importFrom lspline lspline
#' @importFrom latex2exp TeX
#'
#' @examples
#' library(sl3)
#' set.seed(3287)
#' n_obs <- 1000
#' W <- as.numeric(replicate(1, rbinom(n_obs, 1, 0.5)))
#' A <- as.numeric(rnorm(n_obs, mean = 2 * W, sd = 1))
#' Y <- rbinom(n_obs, 1, plogis(2 * A - W))
#' msm <- msm_vimshift(
#'   W = W, A = A, Y = Y, estimator = "tmle",
#'   g_fit_args = list(
#'     fit_type = "sl",
#'     sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
#'   ),
#'   Q_fit_args = list(
#'     fit_type = "glm",
#'     glm_formula = "Y ~ ."
#'   ),
#'   delta_grid = seq(-1, 1, 0.25)
#' )
#' plot(msm)
#'
#' # fit a linear spline with knot at 0
#' set.seed(8293)
#' n_obs <- 1000
#' W <- as.numeric(replicate(1, rbinom(n_obs, 1, 0.5)))
#' A <- as.numeric(rnorm(n_obs, mean = 2 * W, sd = 1))
#' Y <- rbinom(n_obs, 1, plogis(0.1 * A * (A >= 0) - 3 * A * (A < 0) - W))
#' msm <- msm_vimshift(
#'   W = W, A = A, Y = Y, estimator = "tmle",
#'   g_fit_args = list(
#'     fit_type = "sl",
#'     sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
#'   ),
#'   Q_fit_args = list(
#'     fit_type = "glm",
#'     glm_formula = "Y ~ ."
#'   ),
#'   delta_grid = seq(-1, 1, 0.25),
#'   msm_form = list(type = "piecewise", knot = 0)
#' )
#' plot(msm)
#' @export
plot.txshift_msm <- function(x, ...) {
  # create and rename data in MSM fitting for ggplot compatibility
  msm_data <- data.table::copy(x[["param_est"]])
  data.table::setnames(msm_data, c("psi", "delta"), c("y", "x"))

  # extract MSM specs
  msm_type <- x[["msm_type"]]
  msm_knot <- x[["msm_knot"]]

  if (msm_type != "linear" && !is.na(msm_knot)) {
    # fit working MSM regression line if there's a knot point
    msm_mod <- stats::lm(y ~ lspline::lspline(x, msm_knot, marginal = TRUE),
      data = msm_data
    )

    # build geom for MSM in plot
    geom_msm <- ggplot2::geom_smooth(
      method = "lm",
      formula = stats::formula(msm_mod),
      se = FALSE,
      color = "black",
      size = 0.5
    )
  } else {
    # build geom for MSM in plot
    geom_msm <- ggplot2::geom_abline(
      intercept = x[["msm_est"]][["param_est"]][1],
      slope = x[["msm_est"]][["param_est"]][2],
      size = 0.5
    )
  }

  # create plot
  p_msm <- ggplot2::ggplot(data = msm_data, ggplot2::aes_string("x", "y")) +
    ggplot2::geom_point(size = 5) +
    ggplot2::geom_errorbar(
      ggplot2::aes_string(
        ymin = "ci_lwr",
        ymax = "ci_upr"
      ),
      position = "dodge", linetype = "dashed",
      width = 0.05
    ) +
    geom_msm +
    ggplot2::labs(
      x = latex2exp::TeX("Shift in treatment $\\delta$"),
      y = latex2exp::TeX("Counterfactual mean $EY_{A + \\delta(W)}$"),
      title = "Estimated mean counterfactual outcome under shifted treatment",
      subtitle = paste(
        "with marginal confidence intervals and",
        msm_type, "working MSM for summarization"
      )
    ) +
    ggplot2::theme_bw()

  # output plot
  return(p_msm)
}
