#' Iterative Procedure to Compute an IPCW-TMLE
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
#' @param Qn_estim A \code{data.table} corresponding to the outcome regression.
#'  This is produced by invoking the internal function \code{est_Q}.
#' @param Hn_estim A \code{data.table} corresponding to values produced in the
#'  computation of the auxiliary ("clever") covariate. This is produced easily
#'  by invoking the internal function \code{est_Hn}.
#' @param ipcw_mech A \code{numeric} vector containing values that describe the
#'  censoring mechanism for all of the observations. Note that such values are
#'  estimated by regressing the censoring covariates \code{V} on the observed
#'  censoring \code{C} and thus correspond to predicted probabilities of being
#'  censored for each observation.
#' @param ipc_weights_all A \code{numeric} vector of inverse probability of
#'  censoring weights. These are equivalent to \code{C / ipcw_mech} in any
#'  initial run of this function. Updated values of this vector are provided as
#'  part of the output of this function, which may be used in subsequent calls
#'  that allow convergence to a more efficient estimate.
#' @param fluc_method A \code{character} giving the type of regression to be
#'  used in traversing the fluctuation submodel. The choices are "weighted" and
#'  "standard" -- please consult the literature for details on the differences.
#' @param fit_type A \code{character} providing the type of regression to be fit
#'  in estimating the relationship between the variables composing the censoring
#'  mechanism and the efficient influence function. Choose "glm" for generalized
#'  linear models or "sl" for use of the Super Learner algorithm. If selecting
#'  the latter, the final argument \code{sl_lrnrs} must be provided.
#' @param tol_eif A \code{numeric} providing the largest value to be tolerated
#'  as the mean of the efficient influence function.
#' @param sl_lrnrs A \code{Lrnr_sl} object composed of a collection of learners
#'  provided by the \code{sl3} package. This is constructed externally from this
#'  function and provided as input.
#'
#' @importFrom stats var glm qlogis fitted predict as.formula
#' @importFrom data.table as.data.table set copy
#' @importFrom dplyr "%>%"
#' @importFrom sl3 sl3_Task
#'
#' @keywords internal
#'
#' @export
#'
#' @author Nima Hejazi
#' @author David Benkeser
#
ipcw_tmle_proc <- function(data_in, C, V,
                           Qn_estim, Hn_estim,
                           ipcw_mech, ipc_weights_all,
                           fluc_method = c("standard", "weighted"),
                           fit_type = c("glm", "sl"),
                           tol_eif = 1e-9,
                           sl_lrnrs = NULL) {

  # fit logistic regression to fluctuate along the sub-model with NEW WEIGHTS
  fitted_fluc_mod <- fit_fluc(
    Y = data_in$Y,
    Qn_scaled = Qn_estim,
    Hn = Hn_estim,
    ipc_weights = ipc_weights_all[C == 1],
    method = fluc_method
  )

  # for efficiency, need targeting of the censoring mechnanism estimate
  eps_n <- target_pi(Qn = fitted_fluc_mod$Qn_noshift_star,
                     ipc_weights = ipc_weights_all)

  # compute TMLE and EIF using NEW WEIGHTS and UPDATED SUB-MODEL FLUCTUATION
  tmle_eif_out <- tmle_eif_ipcw(
    fluc_fit_out = fitted_fluc_mod,
    Hn = Hn_estim,
    Y = data_in$Y,
    Delta = C,
    ipc_weights = ipc_weights_all[C == 1],
    tol_eif = tol_eif
  )

  # the efficient influence function equation we're solving looks like
  # pi_mech = missing-ness mechanism weights for ALL observations
  # 0 = (C - pi_mech) * \E(f(eif ~ V, subset = (C = 1)) / pi_mech)
  if (fit_type == "sl" & !is.null(sl_lrnrs)) {
    # organize EIF data as data.table
    eif_data <- data.table::copy(data_in) %>%
      subset(., select = names(V))
    data.table::set(eif_data, j = "eif", value = tmle_eif_out$eif[C == 1])

    # make sl3 task from new data.table
    eif_task <- sl3::sl3_Task$new(
      data = eif_data,
      covariates = names(V),
      outcome = "eif",
      outcome_type = "continuous"
    )

    # train an sl3 Super Learner
    fit_eif_sl <- sl_lrnrs$train(eif_task)

    # create task of full data (censoring variables only) for prediction
    eif_pred_task <- sl3::sl3_Task$new(
      data = data.table::set(
        data.table::as.data.table(V), j = "eif",
        value = rep(NaN, unique(lapply(V, length)))
      ),
      covariates = names(V),
      outcome = "eif",
      outcome_type = "continuous"
    )

    # predict from Super Learner on full data
    eif_pred <- fit_eif_sl$predict(eif_pred_task)
  } else {
    # regression model for relationship between censoring variables and EIF
    eif_mod <- stats::glm(
      stats::as.formula("eif ~ ."),
      data = data.table::as.data.table(list(
        eif = tmle_eif_out$eif[C == 1],
        subset(data_in, select = names(V))
      ))
    )

    # compute expectation by invoking the predict method
    eif_pred <- stats::predict(
      eif_mod,
      newdata = data.table::as.data.table(V)
    ) %>%
      as.numeric()
  }

  # fit logistic regression to fluctuate along the sub-model wrt epsilon
  # QUESTION: DO WE WANT TO USE UPDATED VALUES OF CENSORING MECHANISM?
  ipcw_fluc <- stats::glm(
    stats::as.formula("delta ~ -1 + offset(logit_ipcw) + eif_by_ipcw"),
    data = data.table::as.data.table(
      list(
        delta = C,
        logit_ipcw = stats::qlogis(ipcw_mech),
        eif_by_ipcw = I(eif_pred / ipcw_mech)
      )
    ),
    family = "binomial"
  )

  # now, we can obtain P_n^* from the sub-model fluctuation
  ipcw_fluc_pred <- stats::fitted(ipcw_fluc) %>% as.numeric()

  # this is the mean of the second half of the EIF (for censoring bit...)
  ipcw_eif_out <- mean((C - ipcw_fluc_pred) * (eif_pred / ipcw_fluc_pred))

  # sanity check: score of the logistic regression fluctuation model
  ipc_check <- mean((C - ipcw_fluc_pred) * (eif_pred / ipcw_mech))
  stopifnot(ipc_check < tol_eif)

  # so, now we need weights to feed back into the previous steps
  ipc_weights_all <- C / ipcw_fluc_pred

  # as above, compute TMLE and EIF with NEW WEIGHTS but OLD SUBMODEL FLUCTUATION
  tmle_eif_out <- tmle_eif_ipcw(
    fluc_fit_out = fitted_fluc_mod,
    Hn = Hn_estim,
    Y = data_in$Y,
    Delta = C,
    ipc_weights = ipc_weights_all[C == 1],
    tol_eif = tol_eif
  )

  # need to return output such that we can loop over this function
  return(list(
    ipc_weights = ipc_weights_all,
    tmle_eif = tmle_eif_out,
    ipcw_eif = ipcw_eif_out
  ))
}

