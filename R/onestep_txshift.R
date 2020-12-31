#' One-Step Estimate of Counterfactual Mean of Stochastic Shift Intervention
#'
#' @details Invokes the procedure to construct a one-step estimate of the
#'  counterfactual mean under a modified treatment policy.
#'
#' @param data_internal A \code{data.table} constructed internally by a call to
#'  \code{\link{txshift}}. This contains the data elements needed for computing
#'  the one-step estimator.
#' @param C_samp A \code{numeric} indicator for whether a given observation was
#'  included in the second-stage sample, used to compute an IPC-weighted
#'  one-step estimator in cases where two-stage sampling is performed. Default
#'  assumes no censoring due to sampling.
#' @param V The covariates that are used in determining the sampling procedure
#'  that gives rise to censoring. The default is \code{NULL} and corresponds to
#'  scenarios in which there is no censoring (in which case all values in the
#'  preceding argument \code{C} must be uniquely 1. To specify this, pass in a
#'  NAMED \code{list} identifying variables amongst W, A, Y that are thought to
#'  have played a role in defining the sampling/censoring mechanism (C).
#' @param delta A \code{numeric} value indicating the shift in the treatment to
#'  be used in defining the target parameter. This is defined with respect to
#'  the scale of the treatment (A).
#' @param samp_estim An object providing the value of the censoring mechanism
#'  evaluated across the full data. This object is passed in after being
#'  constructed by a call to the internal function \code{\link{est_samp}}.
#' @param gn_cens_weights TODO: document
#' @param Qn_estim An object providing the value of the outcome evaluated after
#'  imposing a shift in the treatment. This object is passed in after being
#'  constructed by a call to the internal function \code{est_Q}.
#' @param Hn_estim An object providing values of the auxiliary ("clever")
#'  covariate, constructed from the treatment mechanism and required for
#'  targeted minimum loss estimation. This object object should be passed in
#'  after being constructed by a call to the internal function \code{est_Hn}.
#' @param eif_reg_type Whether a flexible nonparametric function ought to be
#'  used in the dimension-reduced nuisance regression of the targeting step for
#'  the censored data case. By default, the method used is a nonparametric
#'  regression based on the Highly Adaptive Lasso (from \pkg{hal9001}). Set
#'  this to \code{"glm"} to instead use a simple linear regression model.
#'  In this step, the efficient influence function (EIF) is regressed against
#'  covariates contributing to the censoring mechanism (i.e., EIF ~ V | C = 1).
#' @param samp_fit_args A \code{list} of arguments, all but one of which are
#'  passed to \code{\link{est_samp}}. For details, consult the documentation
#'  for \code{\link{est_samp}}. The first element (i.e., \code{fit_type}) is
#'  used to determine how this regression is fit: "glm" for generalized linear
#'  model, "sl" for a Super Learner, and "external" for a user-specified input
#'  of the form produced by \code{\link{est_samp}}.
#' @param ipcw_efficiency Whether to invoke an augmentation of the IPCW-TMLE
#'  procedure that performs an iterative process to ensure efficiency of the
#'  resulting estimate. The default is \code{TRUE}; set to \code{FALSE} to use
#'  an IPC-weighted loss rather than the IPC-augmented influence function.
#'
#' @importFrom assertthat assert_that
#' @importFrom data.table as.data.table setnames
#' @importFrom stringr str_detect
#' @importFrom Rdpack reprompt
#'
#' @return S3 object of class \code{txshift} containing the results of the
#'  procedure to compute a one-step estimate of the treatment shift parameter.
onestep_txshift <- function(data_internal,
                            C_samp = rep(1, nrow(data_internal)),
                            V = NULL,
                            delta,
                            samp_estim,
                            gn_cens_weights,
                            Qn_estim,
                            Hn_estim,
                            eif_reg_type = c("hal", "glm"),
                            samp_fit_args,
                            ipcw_efficiency = TRUE) {
  # extract and normalize sampling mechanism weights
  samp_weights <- C_samp / samp_estim
  samp_weights_norm <- samp_weights / sum(samp_weights)

  # invoke efficient IPC-weighted one-step if necessary
  if (!all(C_samp == 1) && ipcw_efficiency) {
    # checks for necessary components for augmented EIF procedure
    assertthat::assert_that(!is.null(samp_estim))
    assertthat::assert_that(!is.null(V))

    # quantities to be updated across targeting iterations
    pi_mech_star <- samp_estim
    Qn_estim_updated <- Qn_estim

    # update by submodel fluctuation, re-compute EIF, and update EIF
    os_ipcw_eif <- ipcw_eif_update(
      data_internal = data_internal,
      C_samp = C_samp,
      V = V,
      ipc_mech = pi_mech_star,
      # NOTE: we update pi_mech_star and samp_weights in this procedure, so
      #       only need to rescale by factor gn_cens_weights each iteration
      ipc_weights = (gn_cens_weights * samp_weights[C_samp == 1]),
      Qn_estim = Qn_estim,
      # NOTE: directly pass in Hn since gn never updated in this procedure
      Hn_estim = Hn_estim,
      estimator = "onestep",
      eif_reg_type = eif_reg_type
    )

    # overwrite and update quantities to be used in the next iteration
    Qn_estim_updated <- data.table::as.data.table(
      list(
        # NOTE: need to re-scale estimated outcomes values within bounds of Y
        scale_to_unit(
          vals = os_ipcw_eif$Qn_estim$noshift
        ),
        scale_to_unit(
          vals = os_ipcw_eif$Qn_estim$upshift
        )
      )
    )
    data.table::setnames(Qn_estim_updated, names(Qn_estim))

    # updated sampling weights and stabilize
    samp_weights <- os_ipcw_eif$ipc_weights
    pi_mech_star <- os_ipcw_eif$pi_mech_star

    # compute updated mean of efficient influence function
    eif_ipcw <- os_ipcw_eif$eif_eval$eif - os_ipcw_eif$ipcw_eif_component
    eif_mean <- mean(eif_ipcw)
    eif_var <- var(eif_ipcw) / length(samp_weights)

    # update the one-step estimator with mean of the 2nd half of augmented EIF
    # NOTE: the 2nd half of the EIF is actually a correction (with a minus sign
    #       in front of it) so the mean ought to be subtracted, not added
    psi_os <- os_ipcw_eif$eif_eval$psi - mean(os_ipcw_eif$ipcw_eif_component)

    # create output object
    txshift_out <- unlist(
      list(
        psi = psi_os,
        var = eif_var,
        eif = list(eif_ipcw),
        os_ipcw_eif$eif_eval,
        iter_res = NULL,
        n_iter = 0,
        estimator = "onestep",
        outcome = list(data_internal$Y),
        delta = delta
      ),
      recursive = FALSE
    )
  } else {
    # compute inefficient one-step estimate and EIF
    os_eif_out <- eif(
      Y = data_internal$Y,
      Qn = Qn_estim,
      Hn = Hn_estim,
      estimator = "onestep",
      C_samp = C_samp,
      ipc_weights = (gn_cens_weights * samp_weights[C_samp == 1])
    )

    # create output object
    txshift_out <- unlist(
      list(
        psi = os_eif_out[["psi"]],
        var = os_eif_out[["var"]],
        eif = list(os_eif_out[["eif"]]),
        iter_res = NULL,
        n_iter = 0,
        estimator = "onestep",
        outcome = list(data_internal$Y),
        delta = delta
      ),
      recursive = FALSE
    )
  }

  # S3-ify and return output object
  class(txshift_out) <- "txshift"
  return(txshift_out)
}
