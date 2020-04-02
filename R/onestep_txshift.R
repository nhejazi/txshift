#' Compute One-Step Estimate of Counterfactual Mean Under Stochastic Shift
#' Intervention
#'
#' @details Invokes the procedure to construct a one-step estimate of the
#'  counterfactual mean under a modified treatment policy.
#'
#' @param data_internal A \code{data.table} constructed internally by a call to
#'  \code{\link{txshift}}. This contains the data elements needed for computing
#'  the one-step estimator.
#' @param C A \code{numeric} indicator for whether a given observation censored
#'  in the two-phase sampling procedure, used to compute an IPC-weighted
#'  one-step estimator in cases where two-stage sampling is performed. Default
#'  assumes no censoring.
#' @param V The covariates that are used in determining the sampling procedure
#'  that gives rise to censoring. The default is \code{NULL} and corresponds to
#'  scenarios in which there is no censoring (in which case all values in the
#'  preceding argument \code{C} must be uniquely 1. To specify this, pass in a
#'  NAMED \code{list} identifying variables amongst W, A, Y that are thought to
#'  have played a role in defining the sampling/censoring mechanism (C).
#' @param delta A \code{numeric} value indicating the shift in the treatment to
#'  be used in defining the target parameter. This is defined with respect to
#'  the scale of the treatment (A).
#' @param ipcw_estim An object providing the value of the censoring mechanism
#'  evaluated across the full data. This object is passed in after being
#'  constructed by a call to the internal function \code{est_ipcw}.
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
#' @param ipcw_fit_args A \code{list} of arguments, all but one of which are
#'  passed to \code{\link{est_ipcw}}. For details, consult the documentation
#'  for \code{\link{est_ipcw}}. The first element (i.e., \code{fit_type}) is
#'  used to determine how this regression is fit: "glm" for generalized linear
#'  model, "sl" for a Super Learner, and "fit_spec" a user-specified input of
#'  the form produced by \code{\link{est_ipcw}}. NOTE THAT this first argument
#'  is not passed to \code{\link{est_ipcw}}.
#' @param ipcw_efficiency Whether to invoke an augmentation of the IPCW-TMLE
#'  procedure that performs an iterative process to ensure efficiency of the
#'  resulting estimate. The default is \code{TRUE}; set to \code{FALSE} to use
#'  an IPC-weighted loss rather than the IPC-augment influence function.
#'
#' @importFrom data.table as.data.table setnames
#' @importFrom stringr str_detect
#' @importFrom Rdpack reprompt
#'
#' @return S3 object of class \code{txshift} containing the results of the
#'  procedure to compute a one-step estimate of the treatment shift parameter.
onestep_txshift <- function(data_internal,
                            C = rep(1, nrow(data_internal)),
                            V = NULL,
                            delta,
                            ipcw_estim,
                            Qn_estim,
                            Hn_estim,
                            eif_reg_type = c("hal", "glm"),
                            ipcw_fit_args,
                            ipcw_efficiency = TRUE) {
  # extract and normalize sampling mechanism weights
  cens_weights <- C / ipcw_estim$pi_mech
  cens_weights_norm <- cens_weights / sum(cens_weights)

  # invoke efficient IPC-weighted one-step if necessary
  if (ipcw_efficiency & !all(C == 1) & !is.null(V) & !is.null(ipcw_estim)) {
    # quantities to be updated across iterations
    pi_mech_star <- ipcw_estim$pi_mech
    Qn_estim_updated <- Qn_estim

    # update sub-model fluctuation, re-compute EIF, and update EIF
    os_ipcw_eif <- ipcw_eif_update(
      data_in = data_internal,
      C = C,
      V = V,
      ipc_mech = pi_mech_star,
      ipc_weights = cens_weights,
      ipc_weights_norm = cens_weights_norm,
      Qn_estim = Qn_estim,
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
    cens_weights <- os_ipcw_eif$ipc_weights
    cens_weights_norm <- os_ipcw_eif$ipc_weights_norm
    pi_mech_star <- os_ipcw_eif$pi_mech_star

    # compute updated mean of efficient influence function
    eif_ipcw <- os_ipcw_eif$eif_eval$eif - os_ipcw_eif$ipcw_eif
    eif_ipcw_var <- var(eif_ipcw) / length(C)

    # update the one-step estimator with mean of the 2nd half of augmented EIF
    # NOTE: the 2nd half of the EIF is actually a correction (with a minus sign
    #       in front of it) so the mean ought to be subtracted, not added
    psi_onestep <- os_ipcw_eif$eif_eval$psi - mean(os_ipcw_eif$ipcw_eif)

    # create output object
    txshift_out <- unlist(
      list(
        psi = psi_onestep,
        var = eif_ipcw_var,
        eif = list(eif_ipcw),
        os_ipcw_eif$eif_eval,
        iter_res = NULL,
        n_iter = 0,
        estimator = "onestep",
        outcome = list(data_internal$Y)
      ),
      recursive = FALSE
    )

    # standard (inefficient) one-step of the shift parameter
  } else {
    # compute one-step estimate and EIF for the treatment shift parameter
    os_eif_out <- eif(
      Y = data_internal$Y,
      Qn = Qn_estim,
      Hn = Hn_estim,
      estimator = "onestep",
      Delta = C,
      ipc_weights = cens_weights[C == 1],
      ipc_weights_norm = cens_weights_norm[C == 1]
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
        outcome = list(data_internal$Y)
      ),
      recursive = FALSE
    )
  }

  # S3-ify and return output object
  class(txshift_out) <- "txshift"
  return(txshift_out)
}
