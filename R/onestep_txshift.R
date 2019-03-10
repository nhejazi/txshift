#' Compute Augmented Inverse Probability Weighted Estimate of Counterfactual
#' Mean Under Shifted Treatment
#'
#' @param data_internal ...
#' @param C A \code{numeric} binary vector giving information on whether a given
#'  observation was subject to censoring. This is used to compute an IPCW-AIPW
#'  estimator in cases where two-stage sampling is performed. Default assumes no
#'  censoring was present (i.e., a two-stage design was NOT used). N.B., this is
#'  equivalent to the term %\Delta in the notation used in the original Rose and
#'  van der Laan manuscript that introduced/formulated IPCW-AIPW estimators.
#' @param V The covariates that are used in determining the sampling procedure
#'  that gives rise to censoring. The default is \code{NULL} and corresponds to
#'  scenarios in which there is no censoring (in which case all values in the
#'  preceding argument \code{C} must be uniquely 1. To specify this, pass in a
#'  NAMED \code{list} identifying variables amongst W, A, Y that are thought to
#'  have played a role in defining the sampling/censoring mechanism (C).
#' @param delta A \code{numeric} value indicating the shift in the treatment to
#'  be used in defining the target parameter. This is defined with respect to
#'  the scale of the treatment (A).
#' @param cens_weights ...
#' @param cens_weights_norm ...
#' @param ipcw_estim ...
#' @param Qn_estim ...
#' @param Hn_estim ...
#' @param eif_tol A \code{numeric} giving the convergence criterion for the AIPW
#'  estimator. This is the the maximum mean of the efficient influence function
#'  (EIF) to be used in declaring convergence (theoretically, should be zero).
#' @param max_iter A \code{numeric} integer giving the maximum number of steps
#'  to be taken in iterating to a solution of the efficient influence function.
#' @param eif_reg_type Whether a flexible nonparametric function ought to be
#'  used in the dimension-reduced nuisance regression of the targeting step for
#'  the censored data case. By default, the method used is a nonparametric
#'  regression based on the Highly Adaptive Lasso (from package \code{hal9001}).
#'  Set this to \code{"glm"} to instead use a simple linear regression model.
#'  In this step, the efficient influence function (EIF) is regressed against
#'  covariates contributing to the censoring mechanism (i.e., EIF ~ V | C = 1).
#' @param ipcw_fit_args ...
#' @param ipcw_fit_type ...
#' @param ipcw_efficiency Whether to invoke an augmentation of the IPCW-AIPW
#'  estimation procedure that performs an iterative process to ensure efficiency
#'  of the resulting estimate. The default is \code{TRUE}; only set to
#'  \code{FALSE} if possible inefficiency of the IPCW-AIPW is not a concern.
#'
#' @importFrom data.table as.data.table setnames
#' @importFrom stringr str_detect
#' @importFrom dplyr filter select "%>%"
#' @importFrom Rdpack reprompt
#'
#' @return S3 object of class \code{txshift} containing the results of the
#'  procedure to compute a AIPW estimate of the treatment shift parameter.
#'
#' @export
#
onestep_txshift <- function(data_internal,
                            C = rep(1, nrow(data_internal)),
                            V = NULL,
                            delta,
                            cens_weights,
                            cens_weights_norm,
                            ipcw_estim,
                            Qn_estim,
                            Hn_estim,
                            eif_tol = 1 / nrow(data_internal),
                            max_iter = 1e4,
                            eif_reg_type = c("hal", "glm"),
                            ipcw_fit_args,
                            ipcw_fit_type,
                            ipcw_efficiency = TRUE) {
  # start counter
  n_steps <- 0
  # invoke efficient IPCW-AIPW, per Rose & van der Laan (2011), if necessary
  if (ipcw_efficiency & !all(C == 1) & !is.null(V) & !is.null(ipcw_estim)) {
    # Efficient implementation of the IPCW-AIPW estimator
    eif_mean <- Inf
    conv_res <- replicate(3, rep(NA, max_iter))

    # normalize censoring mechanism weights (to be overwritten)
    cens_weights <- C / ipcw_estim$pi_mech
    cens_weights_norm <- cens_weights / sum(cens_weights)

    # quantities to be updated in iterative procedure (to be overwritten)
    pi_mech_star <- ipcw_estim$pi_mech
    Qn_estim_updated <- Qn_estim

    # iterate procedure until convergence conditions are satisfied
    while (abs(eif_mean) > eif_tol & n_steps < max_iter) {
      # iterate counter
      n_steps <- n_steps + 1

      # update sub-model fluctuation, re-compute EIF, and update EIF
      ipcw_aipw_comp <- ipcw_eif_update(
        data_in = data_internal,
        C = C,
        V = V,
        ipc_mech = pi_mech_star,
        ipc_weights = cens_weights,
        ipc_weights_norm = cens_weights_norm,
        Qn_estim = Qn_estim_updated,
        Hn_estim = Hn_estim,
        estimator = "onestep",
        fit_type = ipcw_fit_type,
        eif_tol = eif_tol,
        sl_lrnrs = ipcw_fit_args$sl_lrnrs,
        eif_reg_type = eif_reg_type
      )

      # overwrite and update quantities to be used in the next iteration
      Qn_estim_updated <- data.table::as.data.table(
        list(
          # NOTE: need to re-scale estimated outcomes values within bounds of Y
          scale_to_unit(
            vals = ipcw_aipw_comp$Qn_estim$noshift
          ),
          scale_to_unit(
            vals = ipcw_aipw_comp$Qn_estim$upshift
          )
        )
      )
      data.table::setnames(Qn_estim_updated, names(Qn_estim))
      cens_weights <- ipcw_aipw_comp$ipc_weights
      cens_weights_norm <- ipcw_aipw_comp$ipc_weights_norm
      pi_mech_star <- ipcw_aipw_comp$pi_mech_star

      # compute updated mean of efficient influence function and save
      eif_mean <- mean(ipcw_aipw_comp$eif_eval$eif - ipcw_aipw_comp$ipcw_eif)
      eif_var <- var(ipcw_aipw_comp$eif_eval$eif - ipcw_aipw_comp$ipcw_eif) /
        length(C)
      conv_res[n_steps, ] <- c(ipcw_aipw_comp$eif_eval$psi, eif_var, eif_mean)
    }
    conv_results <- data.table::as.data.table(conv_res)
    data.table::setnames(conv_results, c("psi", "var", "eif_mean"))

    # replace variance in this object with the updated variance if iterative
    if (exists("eif_var")) {
      ipcw_aipw_comp$eif_eval$var <- eif_var
    }

    # return only the useful convergence results
    conv_results_out <- conv_results[!is.na(rowSums(conv_results)), ]

    # create output object
    txshift_out <- unlist(
      list(
        call = call,
        ipcw_aipw_comp$eif_eval,
        iter_res = list(conv_results_out),
        n_iter = n_steps,
        estimator = "onestep",
        outcome = list(Y)
      ),
      recursive = FALSE
    )

  # standard one-step of the shift parameter / inefficient IPCW-AIPW estimator
  } else {
    # compute one-step estimate and EIF for the treatment shift parameter
    aipw_eif_out <- eif(
      Y = data_internal$Y,
      Qn = Qn_estim,
      Hn = Hn_estim,
      estimator = "onestep",
      Delta = C,
      ipc_weights = cens_weights,
      ipc_weights_norm = cens_weights_norm,
      eif_tol = eif_tol
    )

    # create output object
    txshift_out <- unlist(
      list(
        call = call,
        aipw_eif_out,
        n_iter = n_steps,
        estimator = "onestep",
        outcome = list(Y)
      ),
      recursive = FALSE
    )
  }

  # S3-ify and return output object
  class(txshift_out) <- "txshift"
  return(txshift_out)
}
