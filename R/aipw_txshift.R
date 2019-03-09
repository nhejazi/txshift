#' Compute Augmented Inverse Probability Weighted Estimate of Counterfactual
#' Mean Under Shifted Treatment
#'
#' @param W A \code{matrix}, \code{data.frame}, or similar corresponding to a
#'  set of baseline covariates.
#' @param A A \code{numeric} vector corresponding to a treatment variable. The
#'  parameter of interest is defined as a location shift of this quantity.
#' @param Y A \code{numeric} vector corresponding to an outcome variable.
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
#' @param delta A \code{numeric} value indicating the shift in the treatment to
#'  be used in defining the target parameter. This is defined with respect to
#'  the scale of the treatment (A).
#' @param fluc_method The method to be used in submodel fluctuation step of
#'  the TMLE computation. The choices are "standard" and "weighted".
#' @param eif_tol A \code{numeric} giving the convergence criterion for the TML
#'  estimator. This is the the maximum mean of the efficient influence function
#'  (EIF) to be used in declaring convergence (theoretically, should be zero).
#' @param max_iter A \code{numeric} integer giving the maximum number of steps
#'  to be taken in iterating to a solution of the efficient influence function.
#' @param ipcw_fit_args A \code{list} of arguments, all but one of which are
#'  passed to \code{est_ipcw}. For details, please consult the documentation for
#'  \code{est_ipcw}. The first element of this (i.e., \code{fit_type}) is used
#'  to determine how this regression is fit: "glm" for generalized linear model,
#'  "sl" for a Super Learner, and "fit_spec" a user-specified input of the form
#'  produced by \code{est_ipcw}. NOTE THAT this first argument is not passed to
#'  \code{est_ipcw}.
#' @param g_fit_args A \code{list} of arguments, all but one of which are passed
#'  to \code{est_g}. For further details, please consult the documentation for
#'  \code{est_g}. The first element of this (i.e., \code{fit_type}) is used to
#'  determine how this regression is fit: "glm" for a generalized linear model
#'  for fitting densities via the \code{condensier} package, "sl" for \code{sl3}
#'  learners used to fit Super Learner to densities via \code{Lrnr_condensier},
#'  and "fit_spec" to user-specified input of the form produced by \code{est_g}.
#'  NOTE THAT this first argument is not passed to \code{est_g}.
#' @param Q_fit_args A \code{list} of arguments, all but one of which are passed
#'  to \code{est_Q}. For further details, please consult the documentation for
#'  \code{est_Q}. The first element of this (i.e., \code{fit_type}) is used to
#'  determine how this regression is fit: "glm" for a generalized linear model
#'  for the outcome regression, "sl" for \code{sl3} learners used to fit a Super
#'  Learner for the outcome regression, and "fit_spec" to user-specified input
#'  of the form produced by \code{est_Q}. NOTE THAT this first argument is not
#'  passed to \code{est_g}.
#' @param eif_reg_type Whether a flexible nonparametric function ought to be
#'  used in the dimension-reduced nuisance regression of the targeting step for
#'  the censored data case. By default, the method used is a nonparametric
#'  regression based on the Highly Adaptive Lasso (from package \code{hal9001}).
#'  Set this to \code{"glm"} to instead use a simple linear regression model.
#'  In this step, the efficient influence function (EIF) is regressed against
#'  covariates contributing to the censoring mechanism (i.e., EIF ~ V | C = 1).
#' @param ipcw_efficiency Whether to invoke an augmentation of the IPCW-TMLE
#'  procedure that performs an iterative process to ensure efficiency of the
#'  resulting estimate. The default is \code{TRUE}; only set to \code{FALSE} if
#'  possible inefficiency of the IPCW-TMLE is not a concern.
#' @param ipcw_fit_spec User-specified version of the argument above for fitting
#'  the censoring mechanism (\code{ipcw_fit_args}). Consult the documentation
#'  for that argument for details on how to properly use this. In general, this
#'  should only be used by advanced users familiar with both the underlying
#'  theory and this software implementation of said theory.
#' @param gn_fit_spec User-specified version of the argument above for fitting
#'  the censoring mechanism (\code{g_fit_args}). Consult the documentation for
#'  that argument for details on how to properly use this. In general, this
#'  should only be used by advanced users familiar with both the underlying
#'  theory and this software implementation of said theory.
#' @param Qn_fit_spec User-specified version of the argument above for fitting
#'  the censoring mechanism (\code{Q_fit_args}). Consult the documentation for
#'  that argument for details on how to properly use this. In general, this
#'  should only be used by advanced users familiar with both the underlying
#'  theory and this software implementation of said theory.
#'
#' @importFrom condensier speedglmR6
#' @importFrom data.table as.data.table setnames
#' @importFrom stringr str_detect
#' @importFrom dplyr filter select "%>%"
#' @importFrom Rdpack reprompt
#'
#' @return S3 object of class \code{txshift} containing the results of the
#'  procedure to compute a TML estimate of the treatment shift parameter.
#'
#' @export
#
tmle_txshift <- function(W,
                         A,
                         Y,
                         C = rep(1, length(Y)),
                         V = NULL,
                         delta = 0,
                         fluc_method = c("standard", "weighted"),
                         eif_tol = 1 / length(Y),
                         max_iter = 1e4,
                         ipcw_fit_args = list(
                           fit_type = c("glm", "sl", "fit_spec"),
                           sl_lrnrs = NULL
                         ),
                         g_fit_args = list(
                           fit_type = c("glm", "sl", "fit_spec"),
                           nbins = 35,
                           bin_method = "dhist",
                           bin_estimator =
                             condensier::speedglmR6$new(),
                           parfit = FALSE,
                           sl_lrnrs_dens = NULL
                         ),
                         Q_fit_args = list(
                           fit_type = c("glm", "sl", "fit_spec"),
                           glm_formula = "Y ~ .",
                           sl_lrnrs = NULL
                         ),
                         eif_reg_type = c("hal", "glm"),
                         ipcw_efficiency = TRUE,
                         ipcw_fit_spec = NULL,
                         gn_fit_spec = NULL,
                         Qn_fit_spec = NULL) {
  ##############################################################################
  # invoke efficient IPCW-TMLE, per Rose & van der Laan (2011), if necessary
  ##############################################################################
  n_steps <- 0
  if (ipcw_efficiency & !all(C == 1) & !is.null(V)) {
    # Efficient implementation of the IPCW-TMLE
    eif_mean <- Inf
    conv_res <- replicate(3, rep(NA, max_iter))

    # normalize censoring mechanism weights (to be overwritten)
    cens_weights <- C / ipcw_out$pi_mech
    cens_weights_norm <- cens_weights / sum(cens_weights)

    # quantities to be updated in iterative procedure (to be overwritten)
    pi_mech_star <- ipcw_out$pi_mech
    Qn_estim_use <- Qn_estim

    # iterate procedure until convergence conditions are satisfied
    while (abs(eif_mean) > eif_tol & n_steps < max_iter) {
      # iterate counter
      n_steps <- n_steps + 1

      # update sub-model fluctuation, re-compute EIF, and update EIF
      ipcw_tmle_comp <- ipcw_tmle_proc(
        data_in = data_internal,
        C = C,
        V = V_in,
        ipc_mech = pi_mech_star,
        ipc_weights = cens_weights,
        ipc_weights_norm = cens_weights_norm,
        Qn_estim = Qn_estim_use,
        Hn_estim = Hn_estim,
        fluc_method = fluc_method,
        fit_type = ipcw_fit_type,
        eif_tol = eif_tol,
        sl_lrnrs = ipcw_fit_args$sl_lrnrs,
        eif_reg_type = eif_reg_type
      )

      # overwrite and update quantities to be used in the next iteration
      Qn_estim_use <- data.table::as.data.table(
        list(
          # NOTE: need to re-scale estimated outcomes values within bounds of Y
          bound_scaling(
            Y = Y,
            pred_vals = ipcw_tmle_comp$fluc_mod_out$Qn_noshift_star,
            scale_target =
              ipcw_tmle_comp$fluc_mod_out$Qn_noshift_star,
            scale_type = "bound_in_01"
          ),
          bound_scaling(
            Y = Y,
            pred_vals = ipcw_tmle_comp$fluc_mod_out$Qn_shift_star,
            scale_target =
              ipcw_tmle_comp$fluc_mod_out$Qn_shift_star,
            scale_type = "bound_in_01"
          )
        )
      )
      data.table::setnames(Qn_estim_use, names(Qn_estim))
      cens_weights <- ipcw_tmle_comp$ipc_weights
      cens_weights_norm <- ipcw_tmle_comp$ipc_weights_norm
      pi_mech_star <- ipcw_tmle_comp$pi_mech_star

      # compute updated mean of efficient influence function and save
      eif_mean <- mean(ipcw_tmle_comp$tmle_eif$eif - ipcw_tmle_comp$ipcw_eif)
      eif_var <- var(ipcw_tmle_comp$tmle_eif$eif - ipcw_tmle_comp$ipcw_eif) /
        length(C)
      conv_res[n_steps, ] <- c(ipcw_tmle_comp$tmle_eif$psi, eif_var, eif_mean)
    }
    conv_results <- data.table::as.data.table(conv_res)
    data.table::setnames(conv_results, c("psi", "var", "eif_mean"))
    ##############################################################################
    # standard TMLE of the shift parameter / inefficient IPCW-TMLE
    ##############################################################################
  } else {
    # fit logistic regression to fluctuate along the sub-model
    fitted_fluc_mod <- fit_fluc(
      Y = data_internal$Y,
      Qn_scaled = Qn_estim,
      Hn = Hn_estim,
      ipc_weights = cens_weights,
      method = fluc_method
    )
    # compute TML estimate and EIF for the treatment shift parameter
    tmle_eif_out <- tmle_eif(
      fluc_mod_out = fitted_fluc_mod,
      Hn = Hn_estim,
      Y = data_internal$Y,
      Delta = C,
      ipc_weights = cens_weights,
      ipc_weights_norm = cens_weights_norm,
      tol_eif = eif_tol
    )
  }

  # TODO: what should AIPW output to the wrapper function?
}

