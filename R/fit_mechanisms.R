#' Estimate the Treatment Mechanism (Propensity Score)
#'
#' Compute the propensity score (treatment mechanism) for the observed data,
#' including the shift. This returns the propensity score for the observed data
#' (at A_i) and the shift of interest (at A_i - delta).
#'
#' @param A A \code{numeric} vector of observed treatment values.
#' @param W A \code{numeric} matrix of observed baseline covariate values.
#' @param delta A \code{numeric} value identifying a shift in the observed value
#'  of the treatment under which observations are to be evaluated.
#' @param ipc_weights A \code{numeric} vector of observation-level weights, as
#'  produced by the internal procedure to estimate the censoring mechanism
#'  \code{estimate-ipc_weights}.
#' @param fit_type A \code{character} specifying whether to use Super Learner to
#'  fit the density estimation.
#' @param sl_lrnrs_dens An object containing a set of instantiated learners from
#'  the \code{sl3} package, to be used in fitting an ensemble model.
#' @param std_args A \code{list} of arguments to be passed to \code{fit_density}
#'  from the \code{condensier} package when the argument \code{fit_type} is set
#'  to "glm" (not invoking Super Learner).
#'
#' @importFrom data.table as.data.table setnames set copy
#' @importFrom condensier fit_density predict_probability speedglmR6
#' @importFrom sl3 sl3_Task
#'
#' @export
#'
#
est_g <- function(A,
                  W,
                  delta = 0,
                  ipc_weights = rep(1, length(A)),
                  fit_type = c("sl", "glm"),
                  sl_lrnrs_dens = NULL,
                  std_args = list(
                    nbins = 25,
                    bin_method = "dhist",
                    bin_estimator = condensier::speedglmR6$new(),
                    parfit = FALSE
                  )) {
  ##############################################################################
  # make data objects from inputs
  ##############################################################################
  data_in <- data.table::as.data.table(cbind(A, W))
  if (!is.matrix(W)) W <- as.matrix(W)
  data.table::setnames(data_in, c("A", colnames(W)))
  data.table::set(data_in, j = "ipc_weights", value = ipc_weights)

  # need a data set with the treatment stochastically shifted DOWNWARDS A-delta
  data_in_downshifted <- data.table::copy(data_in)
  data.table::set(data_in_downshifted, j = "A", value = tx_shift(
    A = data_in$A, delta = -delta
  ))

  # need a data set with the treatment stochastically shifted UPWARDS A+delta
  data_in_upshifted <- data.table::copy(data_in)
  data.table::set(data_in_upshifted, j = "A", value = tx_shift(
    A = data_in$A, delta = delta
  ))

  # need a data set with the treatment stochastically shifted UPWARDS A+2delta
  data_in_upupshifted <- data.table::copy(data_in)
  data.table::set(data_in_upupshifted, j = "A", value = tx_shift(
    A = data_in$A, delta = 2 * delta
  ))

  ##############################################################################
  # if fitting sl3 density make sl3 tasks from the data
  ##############################################################################
  if (fit_type == "sl" & !is.null(sl_lrnrs_dens)) {
    # sl3 task for data with treatment UNSHIFTED
    sl_task <- sl3::sl3_Task$new(
      data = data_in,
      outcome = "A",
      covariates = colnames(W),
      weights = "ipc_weights"
    )

    # sl3 task for data with treatment shifted DOWNWARDS A-delta
    sl_task_downshifted <- sl3::sl3_Task$new(
      data = data_in_downshifted,
      outcome = "A",
      covariates = colnames(W),
      weights = "ipc_weights"
    )

    # sl3 task for data with treatment shifted UPWARDS A+delta
    sl_task_upshifted <- sl3::sl3_Task$new(
      data = data_in_upshifted,
      outcome = "A",
      covariates = colnames(W),
      weights = "ipc_weights"
    )

    # sl3 task for data with treatment shifted UPWARDS A+2delta
    sl_task_upupshifted <- sl3::sl3_Task$new(
      data = data_in_upupshifted,
      outcome = "A",
      covariates = colnames(W),
      weights = "ipc_weights"
    )
  }

  ##############################################################################
  # fit conditional densities with condensier
  ##############################################################################
  if (fit_type == "glm" & is.null(sl_lrnrs_dens)) {
    fit_args <- unlist(
      list(
        X = list(colnames(W)),
        Y = "A", weights = "ipc_weights", std_args
      ),
      recursive = FALSE
    )
    fit_args$input_data <- data_in
    fit_g_dens_glm <- do.call(condensier::fit_density, fit_args)
  } else if (fit_type == "sl" & !is.null(sl_lrnrs_dens)) {
    suppressMessages(
      fit_g_dens_sl <- sl_lrnrs_dens$train(sl_task)
    )
  }

  ##############################################################################
  # predict probabilities for the UNSHIFTED data (A = a)
  ##############################################################################
  if (fit_type == "glm" & is.null(sl_lrnrs_dens)) {
    pred_g_A_noshift <-
      condensier::predict_probability(
        model_fit = fit_g_dens_glm,
        newdata = data_in
      )
  } else if (fit_type == "sl" & !is.null(sl_lrnrs_dens)) {
    suppressMessages(
      pred_g_A_noshift <- fit_g_dens_sl$predict()
    )
  }

  ##############################################################################
  # predict probabilities for the DOWNSHIFTED data (A = a - delta)
  ##############################################################################
  if (fit_type == "glm" & is.null(sl_lrnrs_dens)) {
    pred_g_A_downshifted <-
      condensier::predict_probability(
        model_fit = fit_g_dens_glm,
        newdata = data_in_downshifted
      )
  } else if (fit_type == "sl" & !is.null(sl_lrnrs_dens)) {
    suppressMessages(
      pred_g_A_downshifted <- fit_g_dens_sl$predict(sl_task_downshifted)
    )
  }

  ##############################################################################
  # predict probabilities for the UPSHIFTED data (A = a + delta)
  ##############################################################################
  if (fit_type == "glm" & is.null(sl_lrnrs_dens)) {
    pred_g_A_upshifted <-
      condensier::predict_probability(
        model_fit = fit_g_dens_glm,
        newdata = data_in_upshifted
      )
  } else if (fit_type == "sl" & !is.null(sl_lrnrs_dens)) {
    suppressMessages(
      pred_g_A_upshifted <- fit_g_dens_sl$predict(sl_task_upshifted)
    )
  }

  ##############################################################################
  # predict probabilities for the UPSHIFTED data (A = a + 2*delta)
  ##############################################################################
  if (fit_type == "glm" & is.null(sl_lrnrs_dens)) {
    pred_g_A_upupshifted <-
      condensier::predict_probability(
        model_fit = fit_g_dens_glm,
        newdata = data_in_upupshifted
      )
  } else if (fit_type == "sl" & !is.null(sl_lrnrs_dens)) {
    suppressMessages(
      pred_g_A_upupshifted <- fit_g_dens_sl$predict(sl_task_upupshifted)
    )
  }

  ##############################################################################
  # create output data.tables
  ##############################################################################
  out <- data.table::as.data.table(cbind(
    pred_g_A_downshifted,
    pred_g_A_noshift,
    pred_g_A_upshifted,
    pred_g_A_upupshifted
  ))
  data.table::setnames(out, c("downshift", "noshift", "upshift", "upupshift"))
  return(out)
}

################################################################################

#' Estimate the Outcome Regression
#'
#' Compute the outcome regression for the observed data, including with the
#' shift imposed by the intervention. This returns the propensity score for the
#' observed data (at A_i) and the shift of interest (at A_i - delta).
#'
#' @param Y A \code{numeric} vector of observed outcomes.
#' @param A A \code{numeric} vector of observed treatment values.
#' @param W A \code{numeric} matrix of observed baseline covariate values.
#' @param delta A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the treatment \code{A}. This is passed directly to the internal
#'  function \code{tx_shift} and is currently limited to additive shifts.
#' @param ipc_weights A \code{numeric} vector of observation-level weights, as
#'  produced by the internal procedure to estimate the censoring mechanism
#'  \code{estimate-ipc_weights}.
#' @param fit_type A \code{character} indicating whether to use GLMs or Super
#'  Learner to fit the outcome regression. If the option "glm" is selected, the
#'  argument \code{glm_formula} must NOT be \code{NULL}, instead containing a
#'  model formula (in the style of \code{stats::glm}) as a \code{character}. If
#'  the option "sl" is selected, the argument \code{sl_lrnrs} must NOT be
#'  \code{NULL}, instead being a pre-defined object of class \code{Lrnr_sl},
#'  specifying learners and a metalearner for the Super Learner fit. Consult the
#'  documentation of the \code{sl3} package for details on Super Learner fits.
#' @param glm_formula A \code{character} corresponding to a \code{formula} to be
#'  used in fitting a generalized linear model via \code{stats::glm}.
#' @param sl_lrnrs An object containing a set of instantiated learners from the
#'  \code{sl3} package, to be used in fitting an ensemble model.
#'
#' @importFrom stats glm as.formula predict
#' @importFrom data.table as.data.table setnames copy set
#' @importFrom stringr str_detect
#' @importFrom dplyr "%>%"
#' @importFrom sl3 sl3_Task
#'
#' @export
#
est_Q <- function(Y,
                  A,
                  W,
                  delta = 0,
                  ipc_weights = rep(1, length(Y)),
                  fit_type = c("sl", "glm"),
                  glm_formula = "Y ~ .",
                  sl_lrnrs = NULL) {
  ##############################################################################
  # make data objects from inputs but using the outcome y_star instead of y
  ##############################################################################
  # scale the outcome for logit transform
  y_star <- scale_to_unit(vals = Y)

  # generate the data objects for fitting the outcome regression
  data_in <- data.table::as.data.table(cbind(y_star, A, W))
  if (!is.matrix(W)) W <- as.matrix(W)
  data.table::setnames(data_in, c("Y", "A", paste0("W", seq_len(ncol(W)))))
  names_W <- colnames(data_in)[stringr::str_detect(colnames(data_in), "W")]
  data.table::set(data_in, j = "ipc_weights", value = ipc_weights)

  # need a data set with the treatment stochastically shifted UPWARDS...
  data_in_shifted <- data.table::copy(data_in)
  data.table::set(data_in_shifted, j = "A", value = tx_shift(
    A = data_in$A, delta = delta
  ))

  ##############################################################################
  # fit a logistic regression and extract the predicted probabilities
  ##############################################################################
  if (fit_type == "glm" & !is.null(glm_formula)) {
    # need to remove IPCW column from input data.table object for GLM fits
    data.table::set(data_in, j = "ipc_weights", value = NULL)
    data.table::set(data_in_shifted, j = "ipc_weights", value = NULL)

    # fit a logistic regression for the (scaled) outcome
    suppressWarnings(
      fit_Qn <- stats::glm(
        formula = stats::as.formula(glm_formula),
        family = "binomial",
        data = data_in,
        weights = ipc_weights
      )
    )

    # predict Qn for the un-shifted data (A = a)
    pred_star_Qn <- stats::predict(
      object = fit_Qn,
      newdata = data_in,
      type = "response"
    ) %>%
      as.numeric()

    # predict Qn for the shifted data (A = a + delta)
    pred_star_Qn_shifted <- stats::predict(
      object = fit_Qn,
      newdata = data_in_shifted,
      type = "response"
    ) %>%
      as.numeric()
  }

  ##############################################################################
  # fit a binary Super Learner and get the predicted probabilities
  ##############################################################################
  if (fit_type == "sl" & !is.null(sl_lrnrs)) {
    # make sl3 task for original data
    task_noshift <- sl3::sl3_Task$new(
      data = data_in,
      covariates = c("A", names_W),
      outcome = "Y",
      outcome_type = "quasibinomial",
      weights = "ipc_weights"
    )

    # make sl3 task for data with the shifted treatment
    task_shifted <- sl3::make_sl3_Task(
      data = data_in_shifted,
      covariates = c("A", names_W),
      outcome = "Y",
      outcome_type = "quasibinomial",
      weights = "ipc_weights"
    )

    # fit new Super Learner to the no-shift data and predict
    sl_fit_noshift <- sl_lrnrs$train(task_noshift)
    pred_star_Qn <- sl_fit_noshift$predict()

    # predict with Super Learner from unshifted data on the shifted data
    pred_star_Qn_shifted <- sl_fit_noshift$predict(task_shifted)
  }

  ##############################################################################
  # clean up predictions and generate output
  ##############################################################################
  # avoid values that are exactly 0 or 1 in the scaled Qn and Qn_shifted
  pred_star_Qn <- bound_precision(vals = pred_star_Qn)
  pred_star_Qn_shifted <- bound_precision(vals = pred_star_Qn_shifted)

  # create output data frame and return result
  out <- data.table::as.data.table(cbind(pred_star_Qn, pred_star_Qn_shifted))
  data.table::setnames(out, c("noshift", "upshift"))
  return(out)
}

################################################################################

#' Estimate Inverse Probability of Censoring Weights
#'
#' @param V A \code{numeric} vector, \code{matrix}, \code{data.frame} or similar
#'  object giving the observed values of the covariates known to potentially
#'  inform the censoring mechanism.
#' @param Delta A \code{numeric} vector giving observed values of the indicator
#'  function corresponding to the censoring mechanism.
#' @param fit_type A \code{character} indicating whether to perform the fit
#'  using GLMs or a Super Learner. If use of Super Learner is desired, then the
#'  argument \code{sl_lrnrs} must be provided.
#' @param sl_lrnrs An \code{R6} object of class \code{Lrnr_sl}, a Super Learner
#'  object created externally using the \code{sl3} package.
#'
#' @importFrom stats glm predict binomial
#' @importFrom data.table as.data.table setnames
#' @importFrom dplyr "%>%"
#' @importFrom sl3 sl3_Task
#'
#' @return A \code{list} containing a \code{numeric} vector corresponding to the
#'  inverse probability of censoring weights that are required for computing an
#'  IPCW-TMLE and \code{numeric} vector of the estimated missingness mechanism.
#'  Formally, the former is nothing more than %\frac{\Delta}{\Pi_n}, where the
#'  term %\Pi_n is simply the predicted probability of belonging to a censoring
#'  class as computed using standard logistic regression.
#'
#' @export
#
est_ipcw <- function(V,
                     Delta,
                     fit_type = c("sl", "glm"),
                     sl_lrnrs = NULL) {
  ##############################################################################
  # make data objects from inputs
  ##############################################################################
  fit_type <- match.arg(fit_type)
  data_in <- data.table::as.data.table(cbind(Delta, V))
  if (!is.matrix(V)) V <- as.matrix(V)
  names_V <- paste0("V", seq_len(ncol(V)))
  data.table::setnames(data_in, c("Delta", names_V))

  ##############################################################################
  # make sl3 tasks from the data if fitting Super Learners
  ##############################################################################
  if (fit_type == "sl" & !is.null(sl_lrnrs)) {
    sl_task <- sl3::sl3_Task$new(
      data = data_in,
      covariates = names_V,
      outcome = "Delta",
      outcome_type = "binomial"
    )
  }

  ##############################################################################
  # fit logistic regression or binomial SL to get class probabilities for IPCW
  ##############################################################################
  if (fit_type == "glm" & is.null(sl_lrnrs)) {
    ipcw_reg <- stats::glm(
      as.formula("Delta ~ ."),
      family = stats::binomial(),
      data = data_in
    )
    ipcw_probs <- stats::predict(
      object = ipcw_reg,
      newdata = data_in,
      type = "response"
    ) %>%
      as.numeric()
  } else if (fit_type == "sl" & !is.null(sl_lrnrs)) {
    sl_fit <- sl_lrnrs$train(sl_task)
    ipcw_probs <- sl_fit$predict() %>%
      as.numeric()
  }

  # compute the inverse weights as Delta/Pi_n and return this vector
  ipc_weights <- Delta / ipcw_probs
  ipc_weights_out <- ipc_weights[ipc_weights != 0]
  return(list(pi_mech = ipcw_probs, ipc_weights = ipc_weights_out))
}

################################################################################

#' Estimate Auxiliary Covariate from Efficient Influence Function
#'
#' @param gn An estimate of the treatment probability (propensity score), using
#'  the output provided by internal function \code{estimate-propensity_score}.
#' @param a A \code{numeric} vector of the observed values of the treatment.
#' @param w A \code{numeric}, \code{matrix}, \code{data.frame} or similar object
#'  that contains the observed values of the baseline covariates.
#'
#' @importFrom data.table as.data.table setnames
#'
#' @export
#
est_Hn <- function(gn, a = NULL, w = NULL) {
  # set any g(a|w) = 0 values to a very small value above zero
  gn$noshift <- bound_precision(gn$noshift)

  # compute the ratio of the propensity scores for Hn(A,W)
  ratio_g_noshift <- (gn$downshift / gn$noshift) + as.numeric(gn$upshift == 0)

  # compute the ratio of the propensity scores for Hn(d(A,W),W)
  ratio_g_shift <- (gn$noshift / gn$upshift) * as.numeric(gn$upshift != 0) +
    as.numeric(gn$upupshift == 0)

  # set up output table of auxiliary covariate under different shifts
  H_n <- data.table::as.data.table(cbind(ratio_g_noshift, ratio_g_shift))
  data.table::setnames(H_n, c("noshift", "shift"))
  return(H_n)
}

