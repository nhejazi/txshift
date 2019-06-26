#' Compute Targeted Minimum Loss-Based Estimate of Counterfactual Mean Under
#' Stochastically Shifted Treatment
#'
#' This is the primary user-facing wrapper function to be used in invoking the
#' procedure to obtain targeted maximum likelihood / targeted minimum loss-based
#' estimates (TMLE) of the causal parameter discussed in Díaz and van der Laan
#' (2013) and Díaz and van der Laan (2018).
#'
#' @param data_internal A \code{data.table} constructed internally by a call to
#'  \code{txshift}. This contains most of the data elements needed for computing
#'  the targeted maximum likelihood estimator.
#' @param C A \code{numeric} binary vector giving information on whether a given
#'  observation was subject to censoring, used to compute an IPCW-TMLE in cases
#'  where two-stage sampling is performed. Default assumes no censoring.
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
#'  targeted minimum loss-based estimation. This object object should be passed
#'  in after being constructed by a call to the internal function \code{est_Hn}.
#' @param fluc_method The method to be used in submodel fluctuation step of
#'  the TMLE computation. The choices are "standard" and "weighted".
#' @param eif_tol A \code{numeric} giving the convergence criterion for the TML
#'  estimator. This is the maximum mean of the efficient influence function to
#'  be used in declaring convergence.
#' @param max_iter A \code{numeric} integer giving the maximum number of steps
#'  to be taken in iterating to a solution of the efficient influence function.
#' @param eif_reg_type Whether a flexible nonparametric function ought to be
#'  used in the dimension-reduced nuisance regression of the targeting step for
#'  the censored data case. By default, the method used is a nonparametric
#'  regression based on the Highly Adaptive Lasso (from package \code{hal9001}).
#'  Set this to \code{"glm"} to instead use a simple linear regression model.
#'  In this step, the efficient influence function (EIF) is regressed against
#'  covariates contributing to the censoring mechanism (i.e., EIF ~ V | C = 1).
#' @param ipcw_fit_args A \code{list} of arguments, all but one of which are
#'  passed to \code{est_ipcw}. For details, please consult the documentation for
#'  \code{est_ipcw}. The first element of this (i.e., \code{fit_type}) is used
#'  to determine how this regression is fit: "glm" for generalized linear model,
#'  "sl" for a Super Learner, and "fit_spec" a user-specified input of the form
#'  produced by \code{est_ipcw}. NOTE THAT this first argument is not passed to
#'  \code{est_ipcw}.
#' @param ipcw_efficiency Whether to invoke an augmentation of the IPCW-TMLE
#'  procedure that performs an iterative process to ensure efficiency of the
#'  resulting estimate. The default is \code{TRUE}; set to \code{FALSE} to use
#'  an IPC-weighted loss rather than the IPC-augment influence function.
#'
#' @importFrom data.table as.data.table setnames
#' @importFrom stringr str_detect
#' @importFrom dplyr filter select "%>%"
#' @importFrom Rdpack reprompt
#'
#' @return S3 object of class \code{txshift} containing the results of the
#'  procedure to compute a TML estimate of the treatment shift parameter.
#'
#' @export
tmle_txshift <- function(data_internal,
                         C = rep(1, nrow(data_internal)),
                         V = NULL,
                         delta,
                         ipcw_estim,
                         Qn_estim,
                         Hn_estim,
                         fluc_method = c("standard", "weighted"),
                         eif_tol = 1 / nrow(data_internal),
                         max_iter = 1e4,
                         eif_reg_type = c("hal", "glm"),
                         ipcw_fit_args,
                         ipcw_efficiency = TRUE) {
  # initialize counter
  n_steps <- 0

  # normalize censoring mechanism weights (to be overwritten)
  cens_weights <- C / ipcw_estim$pi_mech
  cens_weights_norm <- cens_weights / sum(cens_weights)

  #browser()
  # invoke efficient IPCW-TMLE, per Rose & van der Laan (2011), if necessary
  if (ipcw_efficiency & !all(C == 1) & !is.null(V) & !is.null(ipcw_estim)) {
    # Efficient implementation of the IPCW-TMLE
    eif_mean <- Inf
    conv_res <- matrix(replicate(3, rep(NA_real_, max_iter)), nrow = max_iter)

    # quantities to be updated in iterative procedure (to be overwritten)
    pi_mech_star <- ipcw_estim$pi_mech
    Qn_estim_updated <- Qn_estim

    # iterate procedure until convergence conditions are satisfied
    while (abs(eif_mean) > eif_tol & n_steps < max_iter) {
      # iterate counter
      n_steps <- n_steps + 1

      # update sub-model fluctuation, re-compute EIF, and update EIF
      ipcw_tmle_comp <- ipcw_eif_update(
        data_in = data_internal,
        C = C,
        V = V,
        ipc_mech = pi_mech_star,
        ipc_weights = cens_weights,
        ipc_weights_norm = cens_weights_norm,
        Qn_estim = Qn_estim_updated,
        Hn_estim = Hn_estim,
        estimator = "tmle",
        fluc_method = fluc_method,
        eif_tol = eif_tol,
        eif_reg_type = eif_reg_type
      )

      # overwrite and update quantities to be used in the next iteration
      Qn_estim_updated <- data.table::as.data.table(
        list(
          # NOTE: need to re-scale estimated outcomes values within bounds of Y
          scale_to_unit(
            vals = ipcw_tmle_comp$fluc_mod_out$Qn_noshift_star
          ),
          scale_to_unit(
            vals = ipcw_tmle_comp$fluc_mod_out$Qn_shift_star
          )
        )
      )
      data.table::setnames(Qn_estim_updated, names(Qn_estim))
      cens_weights <- ipcw_tmle_comp$ipc_weights
      cens_weights_norm <- ipcw_tmle_comp$ipc_weights_norm
      pi_mech_star <- ipcw_tmle_comp$pi_mech_star

      # compute updated mean of efficient influence function and save
      eif_mean <- mean(ipcw_tmle_comp$eif_eval$eif - ipcw_tmle_comp$ipcw_eif)
      eif_var <- var(ipcw_tmle_comp$eif_eval$eif - ipcw_tmle_comp$ipcw_eif) /
        length(C)
      conv_res[n_steps, ] <- c(ipcw_tmle_comp$eif_eval$psi, eif_var, eif_mean)
    }
    conv_res <- data.table::as.data.table(conv_res)
    data.table::setnames(conv_res, c("psi", "var", "eif_mean"))

    # replace variance in this object with the updated variance if iterative
    if (exists("eif_var")) {
      ipcw_tmle_comp$eif_eval$var <- eif_var
    }

    # return only the useful convergence results
    conv_res <- conv_res[!is.na(rowSums(conv_res)), ]

    # create output object
    txshift_out <- unlist(
      list(
        call = call,
        ipcw_tmle_comp$eif_eval,
        iter_res = list(conv_res),
        n_iter = n_steps,
        estimator = "tmle",
        outcome = list(data_internal$Y)
      ),
      recursive = FALSE
    )

  # standard TMLE of the shift parameter / inefficient IPCW-TMLE
  } else {
    # fit logistic regression to fluctuate along the sub-model
    fitted_fluc_mod <- fit_fluc(
      Y = data_internal$Y,
      Qn_scaled = Qn_estim,
      Hn = Hn_estim,
      ipc_weights = cens_weights[C == 1],
      method = fluc_method
    )
    # compute TML estimate and EIF for the treatment shift parameter
    tmle_eif_out <- eif(
      Y = data_internal$Y,
      Qn = Qn_estim,
      Hn = Hn_estim,
      estimator = "tmle",
      fluc_mod_out = fitted_fluc_mod,
      Delta = C,
      ipc_weights = cens_weights[C == 1],
      ipc_weights_norm = cens_weights_norm[C == 1],
      eif_tol = eif_tol
    )

    # create output object
    txshift_out <- unlist(
      list(
        call = call,
        tmle_eif_out,
        n_iter = n_steps,
        estimator = "tmle",
        outcome = list(data_internal$Y)
      ),
      recursive = FALSE
    )
  }

  # S3-ify and return output object
  class(txshift_out) <- "txshift"
  return(txshift_out)
}
