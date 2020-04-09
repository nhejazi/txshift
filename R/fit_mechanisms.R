#' Estimate the Generalized Propensity Score (Treatment Mechanism)
#'
#' @details Compute the propensity score (treatment mechanism) for the observed
#'  data, including the shift. This gives the propensity score for the observed
#'  data (at the observed A) and the shift (at A - delta).
#'
#' @param A A \code{numeric} vector of observed treatment values.
#' @param W A \code{numeric} matrix of observed baseline covariate values.
#' @param delta A \code{numeric} value identifying a shift in the observed
#'  value of the treatment under which observations are to be evaluated.
#' @param ipc_weights A \code{numeric} vector of observation-level weights, as
#'  produced by the internal procedure to estimate the censoring mechanism
#'  \code{estimate-ipc_weights}.
#' @param fit_type A \code{character} specifying whether to use Super Learner
#'  (from \pkg{sl3}) or the Highly Adaptive Lasso (from \pkg{hal9001}) to
#'  estimate the conditional treatment density.
#' @param sl_learners_density Object containing a set of instantiated learners
#'  from \pkg{sl3}, to be used in fitting an ensemble model.
#' @param haldensify_args A \code{list} of argument to be directly passed to
#'  \code{\link[haldensify]{haldensify}} when \code{fit_type} is set to
#'  \code{"hal"}. Note that this invokes the Highly Adaptive Lasso instead of
#'  Super Learner and is thus only feasible for relatively small data sets.
#'
#' @importFrom data.table as.data.table setnames set copy
#' @importFrom sl3 sl3_Task
#' @importFrom stats predict
#' @importFrom haldensify haldensify
#' @importFrom assertthat assert_that
#'
#' @return A \code{data.table} with four columns, containing estimates of the
#'  generalized propensity score at a downshift (g(A - delta | W)), no shift
#'  (g(A | W)), an upshift (g(A + delta) | W), and an upshift of magitudie two
#'  (g(A + 2 delta) | W).
est_g <- function(A,
                  W,
                  delta = 0,
                  ipc_weights = rep(1, length(A)),
                  fit_type = c("sl", "hal"),
                  sl_learners_density = NULL,
                  haldensify_args = list(
                    n_bins = c(5, 10),
                    grid_type = c("equal_range", "equal_mass"),
                    lambda_seq = exp(seq(-1, -13, length = 300)),
                    use_future = FALSE
                  )) {
  # set defaults and check arguments
  fit_type <- match.arg(fit_type)
  if (fit_type == "sl") {
    assertthat::assert_that(!is.null(sl_learners_density),
      msg = paste(
        "`sl_learners_density` cannot be NULL",
        "when `fit_type = 'sl'`."
      )
    )
  }

  # make data objects from inputs
  data_in <- data.table::as.data.table(cbind(A, W))
  if (!is.matrix(W)) W <- as.matrix(W)
  data.table::setnames(data_in, c("A", colnames(W)))
  data.table::set(data_in, j = "ipc_weights", value = ipc_weights)

  # need a data set with the treatment stochastically shifted DOWNWARDS A-delta
  data_in_downshifted <- data.table::copy(data_in)
  data.table::set(data_in_downshifted, j = "A", value = shift_additive(
    A = data_in$A, delta = -delta
  ))

  # need a data set with the treatment stochastically shifted UPWARDS A+delta
  data_in_upshifted <- data.table::copy(data_in)
  data.table::set(data_in_upshifted, j = "A", value = shift_additive(
    A = data_in$A, delta = delta
  ))

  # need a data set with the treatment stochastically shifted UPWARDS A+2delta
  data_in_upupshifted <- data.table::copy(data_in)
  data.table::set(data_in_upupshifted, j = "A", value = shift_additive(
    A = data_in$A, delta = 2 * delta
  ))

  # if fitting sl3 density make sl3 tasks from the data
  if (fit_type == "sl" & !is.null(sl_learners_density)) {
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

  # fit conditional densities with haldensify
  if (fit_type == "hal" & is.null(sl_learners_density)) {
    fit_args <- unlist(
      list(
        A = list(A), W = list(W),
        wts = list(ipc_weights),
        haldensify_args
      ),
      recursive = FALSE
    )
    fit_g_dens_hal <- do.call(haldensify::haldensify, fit_args)
  } else if (fit_type == "sl" & !is.null(sl_learners_density)) {
    suppressMessages(
      fit_g_dens_sl <- sl_learners_density$train(sl_task)
    )
  }

  # predict probabilities for the UNSHIFTED data (A = a)
  if (fit_type == "hal" & is.null(sl_learners_density)) {
    pred_g_A_noshift <-
      stats::predict(
        object = fit_g_dens_hal,
        new_A = as.numeric(data_in$A),
        new_W = as.matrix(data_in[, colnames(W), with = FALSE])
      )
  } else if (fit_type == "sl" & !is.null(sl_learners_density)) {
    suppressMessages(
      pred_g_A_noshift <- fit_g_dens_sl$predict()
    )
  }

  # predict probabilities for the DOWNSHIFTED data (A = a - delta)
  if (fit_type == "hal" & is.null(sl_learners_density)) {
    pred_g_A_downshifted <-
      stats::predict(
        object = fit_g_dens_hal,
        new_A = as.numeric(data_in_downshifted$A),
        new_W = as.matrix(data_in[, colnames(W), with = FALSE])
      )
  } else if (fit_type == "sl" & !is.null(sl_learners_density)) {
    suppressMessages(
      pred_g_A_downshifted <- fit_g_dens_sl$predict(sl_task_downshifted)
    )
  }

  # predict probabilities for the UPSHIFTED data (A = a + delta)
  if (fit_type == "hal" & is.null(sl_learners_density)) {
    pred_g_A_upshifted <-
      stats::predict(
        object = fit_g_dens_hal,
        new_A = as.numeric(data_in_upshifted$A),
        new_W = as.matrix(data_in[, colnames(W), with = FALSE])
      )
  } else if (fit_type == "sl" & !is.null(sl_learners_density)) {
    suppressMessages(
      pred_g_A_upshifted <- fit_g_dens_sl$predict(sl_task_upshifted)
    )
  }

  # predict probabilities for the UPSHIFTED data (A = a + 2*delta)
  if (fit_type == "hal" & is.null(sl_learners_density)) {
    pred_g_A_upupshifted <-
      stats::predict(
        object = fit_g_dens_hal,
        new_A = as.numeric(data_in_upupshifted$A),
        new_W = as.matrix(data_in[, colnames(W), with = FALSE])
      )
  } else if (fit_type == "sl" & !is.null(sl_learners_density)) {
    suppressMessages(
      pred_g_A_upupshifted <- fit_g_dens_sl$predict(sl_task_upupshifted)
    )
  }

  # create output data.tables
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
#' @details Compute the outcome regression for the observed data, including
#'  with the shift imposed by the intervention. This returns the propensity
#'  score for the observed data (at A_i) and the shift (at A_i - delta).
#'
#' @param Y A \code{numeric} vector of observed outcomes.
#' @param A A \code{numeric} vector of observed treatment values.
#' @param W A \code{numeric} matrix of observed baseline covariate values.
#' @param delta A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the treatment \code{A}. This is passed to the internal
#'  \code{\link{shift_additive}} and is currently limited to additive shifts.
#' @param ipc_weights A \code{numeric} vector of observation-level weights, as
#'  produced by the internal procedure to estimate the censoring mechanism.
#' @param fit_type A \code{character} indicating whether to use GLMs or Super
#'  Learner to fit the outcome regression. If the option "glm" is selected, the
#'  argument \code{glm_formula} must NOT be \code{NULL}, instead containing a
#'  model formula (as per \code{\link[stats]{glm}}) as a \code{character}. If
#'  the option "sl" is selected, the argument \code{sl_learners} must NOT be
#'  \code{NULL}; instead, an instantiated \code{\link[sl3]{Lrnr_sl}} object,
#'  specifying learners and a metalearner for the Super Learner fit, must be
#'  provided. Consult the documentation of \pkg{sl3} for details.
#' @param glm_formula A \code{character} corresponding to a \code{formula} to
#'  be used in fitting a generalized linear model via \code{\link[stats]{glm}}.
#' @param sl_learners Object containing a set of instantiated learners from the
#'  \pkg{sl3}, to be used in fitting an ensemble model.
#'
#' @importFrom stats glm as.formula predict
#' @importFrom data.table as.data.table setnames copy set
#' @importFrom stringr str_detect
#' @importFrom sl3 sl3_Task
#' @importFrom assertthat assert_that
#'
#' @return A \code{data.table} with two columns, containing estimates of the
#'  outcome mechanism at the natural value of the exposure Q(A, W) and an
#'  upshift of the exposure Q(A + delta, W).
est_Q <- function(Y,
                  A,
                  W,
                  delta = 0,
                  ipc_weights = rep(1, length(Y)),
                  fit_type = c("sl", "glm"),
                  glm_formula = "Y ~ .",
                  sl_learners = NULL) {
  # set defaults and check arguments
  fit_type <- match.arg(fit_type)
  if (fit_type == "sl") {
    assertthat::assert_that(!is.null(sl_learners),
      msg = paste(
        "`sl_learners` cannot be NULL",
        "when `fit_type = 'sl'`."
      )
    )
  }

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
  data.table::set(data_in_shifted, j = "A", value = shift_additive(
    A = data_in$A, delta = delta
  ))

  # fit a logistic regression and extract the predicted probabilities
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
    pred_star_Qn <- unname(stats::predict(
      object = fit_Qn,
      newdata = data_in,
      type = "response"
    ))

    # predict Qn for the shifted data (A = a + delta)
    pred_star_Qn_shifted <- unname(stats::predict(
      object = fit_Qn,
      newdata = data_in_shifted,
      type = "response"
    ))
  }

  # fit a binary Super Learner and get the predicted probabilities
  if (fit_type == "sl" & !is.null(sl_learners)) {
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
    sl_fit_noshift <- sl_learners$train(task_noshift)
    pred_star_Qn <- sl_fit_noshift$predict()

    # predict with Super Learner from unshifted data on the shifted data
    pred_star_Qn_shifted <- sl_fit_noshift$predict(task_shifted)
  }

  # NOTE: clean up predictions and generate output
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
#' @details Compute inverse probability of censoring weights for the two-phase
#'  sampling mechanism. These inverse weights are based on the probability of
#'  appearing in the second-phase sample based on variables measured on all
#'  participants.
#'
#' @param V A \code{numeric} vector, \code{matrix}, \code{data.frame} or
#'  similar object giving the observed values of the covariates known to
#'  potentially inform the censoring mechanism.
#' @param Delta A \code{numeric} vector giving observed values of the indicator
#'  function corresponding to the censoring mechanism.
#' @param fit_type A \code{character} indicating whether to perform the fit
#'  using GLMs or a Super Learner. If use of Super Learner is desired, then the
#'  argument \code{sl_learners} must be provided.
#' @param sl_learners An \code{\link[sl3]{Lrnr_sl}} object, a Super Learner
#'  instantiated externally using \pkg{sl3}.
#'
#' @importFrom stats glm predict binomial
#' @importFrom data.table as.data.table setnames
#' @importFrom sl3 sl3_Task
#' @importFrom assertthat assert_that
#'
#' @return A \code{list} containing a \code{numeric} vector corresponding to
#'  the inverse probability of censoring weights required for computing an
#'  IPCW-TMLE and \code{numeric} vector of the estimated missingness mechanism.
#'  Formally, the former is nothing more than %\frac{\Delta}{\Pi_n}, where the
#'  term %\Pi_n is simply the predicted probability of belonging to a censoring
#'  class as computed using standard logistic regression.
est_ipcw <- function(V,
                     Delta,
                     fit_type = c("sl", "glm"),
                     sl_learners = NULL) {
  # set defaults and check arguments
  fit_type <- match.arg(fit_type)
  if (fit_type == "sl") {
    assertthat::assert_that(!is.null(sl_learners),
      msg = paste(
        "`sl_learners` cannot be NULL",
        "when `fit_type = 'sl'`."
      )
    )
  }

  # make data objects from inputs
  fit_type <- match.arg(fit_type)
  data_in <- data.table::as.data.table(cbind(Delta, V))
  if (!is.matrix(V)) V <- as.matrix(V)
  names_V <- paste0("V", seq_len(ncol(V)))
  data.table::setnames(data_in, c("Delta", names_V))

  # make sl3 tasks from the data if fitting Super Learners
  if (fit_type == "sl" & !is.null(sl_learners)) {
    sl_task <- sl3::sl3_Task$new(
      data = data_in,
      covariates = names_V,
      outcome = "Delta",
      outcome_type = "binomial"
    )
  }

  # fit logistic regression or binomial SL to get class probabilities for IPCW
  if (fit_type == "glm" & is.null(sl_learners)) {
    ipcw_reg <- stats::glm(
      as.formula("Delta ~ ."),
      family = stats::binomial(),
      data = data_in
    )
    ipcw_probs <- unname(stats::predict(
      object = ipcw_reg,
      newdata = data_in,
      type = "response"
    ))
  } else if (fit_type == "sl" & !is.null(sl_learners)) {
    sl_fit <- sl_learners$train(sl_task)
    ipcw_probs <- as.numeric(sl_fit$predict())
  }

  # compute the inverse weights as Delta/Pi_n and return this vector
  ipc_weights <- Delta / ipcw_probs
  ipc_weights_out <- ipc_weights[ipc_weights != 0]
  return(list(pi_mech = ipcw_probs, ipc_weights = ipc_weights_out))
}

################################################################################

#' Estimate Auxiliary Covariate from Efficient Influence Function
#'
#' @details Compute an estimate of the auxiliary covariate required to update
#'  initial estimates via logistic tilting models in targeted minimum loss
#'  estimation.
#'
#' @param gn An estimate of the treatment probability (propensity score), using
#'  the output provided by internal function \code{estimate-propensity_score}.
#'
#' @importFrom data.table as.data.table setnames
#'
#' @return A \code{data.table} with two columns, containing estimates of the
#'  auxiliary covariate at the natural value of the exposure H(A, W) and at the
#'  shifted value of the exposure H(A + delta, W).
est_Hn <- function(gn) {
  # set any g(a|w) = 0 values to a very small value above zero
  gn$noshift <- bound_propensity(gn$noshift)

  # compute the ratio of the propensity scores for Hn(A,W)
  ratio_g_noshift <- (gn$downshift / gn$noshift) + as.numeric(gn$upshift == 0)

  # compute the ratio of the propensity scores for Hn(d(A,W),W)
  ratio_g_shift <- (gn$noshift / gn$upshift) * as.numeric(gn$upshift != 0) +
    as.numeric(gn$upupshift == 0)
  ratio_g_shift[is.na(ratio_g_shift)] <- 1

  # set up output table of auxiliary covariate under different shifts
  H_n <- data.table::as.data.table(cbind(ratio_g_noshift, ratio_g_shift))
  data.table::setnames(H_n, c("noshift", "shift"))
  return(H_n)
}
