#' Test for a trend in the causal effects of a grid of shift interventions
#'
#' @param Y A \code{numeric} vector of the observed outcomes.
#' @param A A \code{numeric} vector corresponding to a treatment variable. The
#'  parameter of interest is defined as a location shift of this quantity.
#' @param W A \code{matrix}, \code{data.frame}, or similar corresponding to a
#'  set of baseline covariates.
#' @param C A \code{numeric} binary vector giving information on whether a given
#'  observation was subject to censoring. This is used to compute an IPCW-TMLE
#'  in cases where two-stage sampling is performed. The default assumes that no
#'  censoring was present (i.e., a two-stage design was NOT used). N.B., this is
#'  equivalent to the term %\Delta in the notation used in the original Rose and
#'  van der Laan manuscript that introduced/formulated IPCW-TML estimators.
#' @param V The covariates that are used in determining the sampling procedure
#'  that gives rise to censoring. The default is \code{NULL} and corresponds to
#'  scenarios in which there is no censoring (in which case all values in the
#'  preceding argument \code{C} must be uniquely 1. To specify this, pass in a
#'  NAMED \code{list} identifying variables amongst W, A, Y that are thought to
#'  have played a role in defining the sampling/censoring mechanism (C).
#' @param delta_grid A \code{numeric} vector giving the individual values of the
#'  shift parameter used in computing each of the TML estimates.
#' @param estimator The type of estimator to be fit, either \code{"tmle"} for
#'  targeted maximum likelihood estimation or \code{"onestep"} for a one-step
#'  augmented inverse probability weighted (AIPW) estimator.
#' @param weighting Whether to weight each parameter estimate by the inverse of
#'  its variance (in order to improve stability of the resultant MSM fit) or to
#'  simply weight all parameter estimates equally.
#' @param ci_level A \code{numeric} indicating the desired coverage level of the
#'  confidence interval to be computed.
#' @param ci_type Whether to construct a simultaneous confidence band covering
#'  all parameter estimates at once or marginal confidence intervals covering
#'  each parameter estimate separately.
#' @param ... Additional arguments to be passed to \code{txshift}.
#'
#' @importFrom stats cov qnorm pnorm
#' @importFrom assertthat assert_that
#' @importFrom tibble as_tibble
#'
#' @export
#
#msm_vimshift <- function(Y,
                         #A,
                         #W,
                         #C,
                         #V,
                         #delta_grid,
                         #estimator = c("tmle", "onestep"),
                         #weighting = c("identity", "variance"),
                         #ci_level = 0.95,
                         #ci_type = c("simultaneous", "marginal"),
                         #...) {

  ## make sure more than one parameter has been estimated for trend
  #assertthat::assert_that(length(tmle_fit_estimates) > 1)

  ## matrix of EIF(O_i) values and estimates across each parameter estimated
  #eif_mat <- sapply(tmle_fit_estimates, `[[`, "eif")
  #psi_vec <- sapply(tmle_fit_estimates, `[[`, "psi")

  ## set weights to be the inverse of the variance of each TML estimate
  #if (is.null(weights)) {
    #weights <- as.numeric(1 / diag(stats::cov(eif_mat)))
  #}

  ## multiplier for CI construction
  #ci_mult <- (c(1, -1) * stats::qnorm((1 - level) / 2))

  ## compute the MSM parameters
  #intercept <- rep(1, length(delta_grid))
  #x_mat <- cbind(intercept, delta_grid)
  #omega <- diag(weights)
  #s_mat <- solve(t(x_mat) %*% omega %*% x_mat) %*% t(x_mat) %*% omega
  #msm_param <- as.vector(s_mat %*% psi_vec)

  ## compute inference for MSM based on individual EIF(O_i) for each parameter
  #msm_eif <- t(tcrossprod(s_mat, eif_mat))
  #msm_var <- diag(stats::cov(msm_eif))
  #msm_se <- sqrt(msm_var / nrow(msm_eif))

  ## build confidence intervals and hypothesis tests for EIF(msm)
  #ci_msm_param <- msm_se %*% t(ci_mult) + msm_param
  #pval_msm_param <- 2 * stats::pnorm(-abs(msm_param / msm_se))

  ## table for output
  #txshift_out <- list(
    #ci_low = psi_vec + ci_mult[1] * sqrt(diag(stats::cov(eif_mat))),
    #psi = psi_vec,
    #ci_high = psi_vec + ci_mult[2] * sqrt(diag(stats::cov(eif_mat)))
  #) %>%
  #tibble::as_tibble()

  #msm_out <- list(
    #param = names(msm_se),
    #ci_low = ci_msm_param[, 1],
    #param_est = msm_param,
    #ci_high = ci_msm_param[, 2],
    #param_se = msm_se,
    #p_value = pval_msm_param
  #) %>%
  #tibble::as_tibble()
  #return(out)
#}
