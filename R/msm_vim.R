#' Marginal structural model for the causal effects of a grid of interventions
#'
#' @details TODO
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
#' @importFrom stats cov qnorm pnorm
#' @importFrom mvtnorm qmvnorm
#'
#' @return TODO
#'
#' @examples
#' n_obs <- 50
#' W <- as.numeric(replicate(1, rbinom(n_obs, 1, 0.5)))
#' A <- as.numeric(rnorm(n_obs, mean = 2 * W, sd = 1))
#' Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
#' msm <- msm_vimshift(
#'   W = W, A = A, Y = Y, estimator = "tmle",
#'   g_fit_args = list(
#'     fit_type = "hal",
#'     n_bins = 3,
#'     grid_type = "equal_mass",
#'     lambda_seq = exp(seq(-1, -9,
#'       length = 50
#'     ))
#'   ),
#'   Q_fit_args = list(
#'     fit_type = "glm",
#'     glm_formula = "Y ~ ."
#'   )
#' )
#' @export
msm_vimshift <- function(Y,
                         A,
                         W,
                         C = rep(1, length(Y)),
                         V = NULL,
                         delta_grid = seq(-0.5, 0.5, 0.5),
                         estimator = c("tmle", "onestep"),
                         weighting = c("identity", "variance"),
                         ci_level = 0.95,
                         ci_type = c("marginal", "simultaneous"),
                         ...) {
  # set default values of arguments
  estimator <- match.arg(estimator)
  weighting <- match.arg(weighting)
  ci_type <- match.arg(ci_type)

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
      return(est)
    })

  # matrix of EIF(O_i) values and estimates across each parameter estimated
  psi_vec <- sapply(est_over_grid, `[[`, "psi")
  eif_mat <- sapply(est_over_grid, `[[`, "eif")

  # set weights to be the inverse of the variance of each TML estimate
  if (weighting == "variance") {
    weights <- as.numeric(1 / diag(stats::cov(eif_mat)))
  } else {
    weights <- rep(1, length(psi_vec))
  }

  # compute the MSM parameters
  intercept <- rep(1, length(delta_grid))
  x_mat <- cbind(intercept, delta_grid)
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
      ci_low = psi_vec + ci_mult[1] * sqrt(diag(stats::cov(eif_mat))),
      psi = psi_vec,
      ci_high = psi_vec + ci_mult[2] * sqrt(diag(stats::cov(eif_mat)))
    )
  )

  # create summary table for MSM estimates
  msm_out <- data.table::as.data.table(list(
    param = names(msm_se),
    ci_low = ci_msm_param[, 1],
    param_est = msm_param,
    ci_high = ci_msm_param[, 2],
    param_se = msm_se,
    p_value = pval_msm_param
  ))

  # complete output for MSM
  out <- list(
    param_est = vimshift_out,
    msm_est = msm_out
  )
  return(out)
}
