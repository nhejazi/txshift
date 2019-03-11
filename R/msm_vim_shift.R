#' Test for a trend in the causal effects of a grid of shift interventions
#'
#' @param Y ...
#' @param A ...
#' @param W ...
#' @param C ...
#' @param delta_grid A \code{numeric} vector giving the individual values of the
#'  shift parameter used in computing each of the TML estimates.
#' @param estimator ...
#' @param weighting A \code{numeric} vector indicating the weights (if any) to
#'  be applied to each of the estimated mean counterfactual outcomes under
#'  posited values of the shift parameter delta. The default is to weight each
#'  estimate by the inverse of its variance, in order to improve stability;
#'  however, this may be changed depending on the exact choice of shift
#'  function.
#' @param ci_level The nominal coverage probability of the confidence interval.
#' @param ci_type ...
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
