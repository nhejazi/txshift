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
#' @importFrom ggplot2 ggplot geom_point geom_errorbar geom_segment geom_smooth
#'  aes aes_string labs theme_bw
#' @importFrom stats formula
#' @importFrom latex2exp TeX
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
#'   plot(msm)
#'
#'   # fit a linear spline with knot at 0
#'   set.seed(8293)
#'   n_obs <- 1000
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
#'   plot(msm)
#' }
#' @export
plot.txshift_msm <- function(x, ...) {
  # build geom for MSM in plot
  if (x$.msm_type == "piecewise") {
    geom_msm <- ggplot2::geom_smooth(
      method = "lm",
      formula = stats::formula(x$msm_fit),
      se = FALSE,
      color = "black",
      size = 0.5,
      linetype = "dashed"
    )
  } else if (x$.msm_type == "linear") {
    intercept <- x$msm_est$param_est[1]
    slope <- x$msm_est$param_est[2]
    delta_grid <- x$.delta_grid
    geom_msm <- ggplot2::geom_segment(
      ggplot2::aes(
        x = min(delta_grid),
        xend = max(delta_grid),
        y = intercept + min(delta_grid) * slope,
        yend = intercept + max(delta_grid) * slope
      ),
      size = 0.5, color = "black", linetype = "dashed"
    )
  }

  # error bars for marginal CIs but band for simultaneous CIs
  if (x$.ci_type == "marginal") {
    geom_ci <- ggplot2::geom_errorbar(
      data = x$.msm_data,
      ggplot2::aes_string(
        ymin = "ci_lwr",
        ymax = "ci_upr"
      ),
      linetype = "dotted",
      width = 0.05
    )
  } else if (x$.ci_type == "simultaneous") {
    geom_ci <- ggplot2::geom_ribbon(
      data = x$.msm_data,
      ggplot2::aes_string(
        ymin = "ci_lwr",
        ymax = "ci_upr"
      ),
      fill = "grey",
      alpha = 0.3
    )
  }

  # create plot
  p_msm <- ggplot2::ggplot(
    data = x$.msm_data,
    ggplot2::aes_string("x", "y")
  ) +
    geom_msm +
    ggplot2::geom_point(size = 3, alpha = 0.75) +
    geom_ci +
    ggplot2::labs(
      x = latex2exp::TeX("Shift in treatment $\\delta$"),
      y = latex2exp::TeX("Counterfactual mean $EY_{A + \\delta(W)}$"),
      title = "Estimated mean counterfactual outcome under shifted treatment",
      subtitle = paste(
        "with", x$.ci_type, "confidence intervals and",
        x$.msm_type, "working MSM for summarization"
      )
    ) +
    ggplot2::theme_bw()

  # output plot
  return(p_msm)
}
