#' Iterative IPCW Update Procedure of Efficient Influence Function
#'
#' An adaptation of the general IPCW-TMLE formulation of Rose & van der Laan as
#' well as its associated algorithm. This may be used to iteratively construct
#' an efficient IPCW-TMLE, which is computed by applying IPC weights to the
#' efficient influence function of the parameter and subsequently updating the
#' IPC weights. The efficient IPCW-TMLE estimate is provided by looping over
#' this function until convergence criteria defined over the efficient influence
#' function are satisfied (e.g., a mean of very nearly zero). For a detailed
#' description of the IPCW-TMLE, please consult Rose and van der Laan (2011)
#' <doi:10.2202/1557-4679.1217>. INTERNAL USE ONLY.
#'
#' @param data_in A \code{data.table} containing variables and observations that
#'  represent the fully observed data. That is, this is a table of the data as
#'  it is seen after a potential censoring process is applied.
#' @param C A \code{numeric} binary vector giving the censoring status of a
#'  given observation.
#' @param V A \code{data.table} giving the values across all observations of all
#'  variables that play a role in the censoring mechanism.
#' @param ipc_mech A \code{numeric} vector containing values that describe the
#'  censoring mechanism for all of the observations. Note that such values are
#'  estimated by regressing the censoring covariates \code{V} on the observed
#'  censoring \code{C} and thus correspond to predicted probabilities of being
#'  censored for each observation.
#' @param ipc_weights A \code{numeric} vector of inverse probability of
#'  censoring weights. These are equivalent to \code{C / ipc_mech} in any
#'  initial run of this function. Updated values of this vector are provided as
#'  part of the output of this function, which may be used in subsequent calls
#'  that allow convergence to a more efficient estimate.
#' @param ipc_weights_norm A \code{numeric} vector of the same weights described
#'  in the previous argument, except normalized in this case.
#' @param Qn_estim A \code{data.table} corresponding to the outcome regression.
#'  This is produced by invoking the internal function \code{est_Q}.
#' @param Hn_estim A \code{data.table} corresponding to values produced in the
#'  computation of the auxiliary ("clever") covariate. This is produced easily
#'  by invoking the internal function \code{est_Hn}.
#' @param estimator The type of estimator to be fit, either \code{"tmle"} for
#'  targeted maximum likelihood estimation or \code{"onestep"} for a one-step
#'  or augmented inverse probability weighted (AIPW) estimator.
#' @param fluc_method A \code{character} giving the type of regression to be
#'  used in traversing the fluctuation submodel. The choices are "weighted" and
#'  "standard" -- please consult the literature for details on the differences.
#' @param eif_tol A \code{numeric} providing the largest value to be tolerated
#'  as the mean of the efficient influence function.
#' @param eif_reg_type Whether a flexible nonparametric function ought to be
#'  used in the dimension-reduced nuisance regression of the targeting step for
#'  the censored data case. By default, the method used is a nonparametric
#'  regression based on the Highly Adaptive Lasso (from package \code{hal9001}).
#'  Set this to \code{"glm"} to instead use a simple linear regression model.
#'  In this step, the efficient influence function (EIF) is regressed against
#'  covariates contributing to the censoring mechanism (i.e., EIF ~ V | C = 1).
#'
#' @importFrom stats var glm qlogis fitted predict as.formula
#' @importFrom data.table as.data.table set copy
#' @importFrom assertthat assert_that
#' @importFrom dplyr "%>%" select
#' @importFrom hal9001 fit_hal
#'
#' @keywords internal
#'
#' @export
#'
#' @author Nima Hejazi
#' @author David Benkeser
ipcw_eif_update <- function(data_in,
                            C,
                            V,
                            ipc_mech,
                            ipc_weights,
                            ipc_weights_norm,
                            Qn_estim,
                            Hn_estim,
                            estimator = c("tmle", "onestep"),
                            fluc_method = NULL,
                            eif_tol = 1e-9,
                            eif_reg_type = c("hal", "glm")) {
  # perform submodel fluctuation if computing TMLE
  if (estimator == "tmle" & !is.null(fluc_method)) {
    # fit logistic regression for submodel fluctuation with updated weights
    fitted_fluc_mod <- fit_fluc(
      Y = data_in$Y,
      Qn_scaled = Qn_estim,
      Hn = Hn_estim,
      ipc_weights = ipc_weights[C == 1],
      method = fluc_method
    )
  } else if (estimator == "onestep" & is.null(fluc_method)) {
    fitted_fluc_mod <- NULL
  }

  # NOTE: the augmented EIF estimating equation being solved takes the form
  # 0 = (C - pi) * \E(f(eif ~ V, subset = (C = 1)) / pi)

  # compute EIF using updated weights and updated fluctuation (if TMLE)
  # NOTE: for one-step, this adds the first half of the EIF as the correction
  #       SO the second half (from the reduced regression) is still needed...
  eif_eval <- eif(
    Y = data_in$Y,
    Qn = Qn_estim,
    Hn = Hn_estim,
    estimator = estimator,
    fluc_mod_out = fitted_fluc_mod,
    Delta = C,
    ipc_weights = ipc_weights[C == 1],
    ipc_weights_norm = ipc_weights_norm[C == 1],
    eif_tol = eif_tol
  )

  # organize EIF data as data.table
  eif_data <- data_in %>%
    data.table::copy() %>%
    dplyr::select(names(V)) %>%
    data.table::as.data.table() %>%
    data.table::set(j = "eif", value = eif_eval$eif[C == 1])

  # estimate the EIF nuisance regression using HAL
  if (eif_reg_type == "hal") {
    # make sl3 task from new data.table
    # eif_task <- sl3::sl3_Task$new(
    # data = eif_data,
    # covariates = names(V),
    # outcome = "eif",
    # outcome_type = "continuous"
    # )

    # train HAL via sl3 learner
    # hal_lrnr <- sl3::Lrnr_hal9001$new()
    # fit_eif_hal <- hal_lrnr$train(eif_task)

    # create task of full data (censoring variables only) for prediction
    # eif_pred_task <- sl3::sl3_Task$new(
    # data = data.table::set(
    # data.table::as.data.table(V),
    # j = "eif",
    # value = rep(0, unique(lapply(V, length)))
    # ),
    # covariates = names(V),
    # outcome = "eif",
    # outcome_type = "continuous"
    # )

    # predict from HAL on full data
    # eif_pred_hal <- fit_eif_hal$predict(eif_pred_task)

    # if flexibility specified, just fit a HAL regression
    X_in <- eif_data %>%
      data.table::copy() %>%
      dplyr::select(names(V)) %>%
      as.matrix()
    colnames(X_in) <- NULL

    # fit HAL with custom arguments to get results similar to halplus
    # NOTE: this produces a warning that could be leading to serious issues
    eif_mod <- hal9001::fit_hal(
      X = X_in,
      Y = as.numeric(eif_data$eif),
      standardize = FALSE,
      fit_type = "glmnet", # "lassi",
      lambda = exp(seq(3, -50, length = 2000)),
      return_lasso = TRUE,
      yolo = FALSE
    )
    # compute expectation by invoking the predict method
    V_pred_in <- as.matrix(V)
    colnames(V_pred_in) <- NULL
    eif_pred <- stats::predict(
      eif_mod,
      new_data = V_pred_in
    ) %>%
      as.numeric()
  } else if (eif_reg_type == "glm") {
    # NOTE: change of formula to examine interactions (for debugging)
    eif_mod <- stats::glm(
      stats::as.formula("eif ~ .^2"),
      data = eif_data
    )
    # compute expectation by invoking the predict method
    eif_pred <- stats::predict(
      eif_mod,
      newdata = V
    ) %>%
      as.numeric()
  }

  # TMLE: fit logistic regression to fluctuate along submodel wrt epsilon
  if (estimator == "tmle") {
    ipcw_fluc_reg_data <-
      data.table::as.data.table(
        list(
          delta = C,
          logit_ipcw = stats::qlogis(bound_precision(ipc_mech)),
          eif_by_ipcw = (eif_pred / bound_precision(ipc_mech))
        )
      )
    # fit fluctuation regression
    ipcw_fluc <- stats::glm(
      stats::as.formula("delta ~ -1 + offset(logit_ipcw) + eif_by_ipcw"),
      data = ipcw_fluc_reg_data,
      family = "binomial"
    )
    # now, we can obtain Pn* from the sub-model fluctuation
    ipcw_pred <- stats::fitted(ipcw_fluc) %>%
      as.numeric()
  } else {
    # just use the initial estimates of censoring probability for one-step
    ipcw_pred <- ipc_mech
  }

  # this is the mean of the second half of the IPCW-EIF
  ipcw_eif_out <- (C - ipcw_pred) * (eif_pred / ipcw_pred)

  # so, now we need weights to feed back into the previous steps
  ipc_weights <- C / ipcw_pred
  ipc_weights_norm <- ipc_weights / sum(ipc_weights)

  # as above, compute TMLE and EIF with NEW WEIGHTS and SUBMODEL FLUCTUATION
  # NOTE: since this is meant to update the EIF components based on the TMLE
  #       update steps, it accomplishes _literally_ nothing for the one-step
  if (estimator == "tmle") {
    eif_eval <- eif(
      Y = data_in$Y,
      Qn = Qn_estim,
      Hn = Hn_estim,
      estimator = estimator,
      fluc_mod_out = fitted_fluc_mod,
      Delta = C,
      ipc_weights = ipc_weights[C == 1],
      ipc_weights_norm = ipc_weights_norm[C == 1],
      eif_tol = eif_tol
    )
  }

  # need to return output such that we can loop over this function
  out <- list(
    Qn_estim = Qn_estim,
    fluc_mod_out = fitted_fluc_mod,
    pi_mech_star = ipcw_pred,
    ipc_weights = ipc_weights,
    ipc_weights_norm = ipc_weights_norm,
    eif_eval = eif_eval,
    ipcw_eif = ipcw_eif_out
  )
  return(out)
}
