#' Estimate the Exposure Mechanism via Generalized Propensity Score
#'
#' @details Compute the propensity score (exposure mechanism) for the observed
#'  data, including the shift. This gives the propensity score for the observed
#'  data (at the observed A) the counterfactual shifted exposure levels (at
#'  {A - delta}, {A + delta}, and {A + 2 * delta}).
#'
#' @param A A \code{numeric} vector of observed exposure values.
#' @param W_g A \code{numeric} matrix of observed baseline covariate values for estimating g.
#' @param delta A \code{numeric} value identifying a shift in the observed
#'  value of the exposure under which observations are to be evaluated.
#' @param samp_weights A \code{numeric} vector of observation-level sampling
#'  weights, as produced by the internal procedure to estimate the two-phase
#'  sampling mechanism \code{\link{est_samp}}.
#' @param fit_type A \code{character} specifying whether to use Super Learner
#'  (from \pkg{sl3}) or the Highly Adaptive Lasso (from \pkg{hal9001}) to
#'  estimate the conditional exposure density.
#' @param sl_learners_density Object containing a set of instantiated learners
#'  from \pkg{sl3}, to be used in fitting an ensemble model.
#' @param haldensify_args A \code{list} of argument to be directly passed to
#'  \code{\link[haldensify]{haldensify}} when \code{fit_type} is set to
#'  \code{"hal"}. Note that this invokes the Highly Adaptive Lasso instead of
#'  Super Learner and is thus only feasible for relatively small data sets.
#'
#' @importFrom data.table as.data.table setnames set copy
#' @importFrom stats predict
#' @importFrom haldensify haldensify
#' @importFrom assertthat assert_that
#'
#' @return A \code{data.table} with four columns, containing estimates of the
#'  generalized propensity score at a downshift (g(A - delta | W)), no shift
#'  (g(A | W)), an upshift (g(A + delta) | W), and an upshift of magnitude two
#'  (g(A + 2 delta) | W).
est_g_exp <- function(A,
                      W_g,
                      delta = 0,
                      samp_weights = rep(1, length(A)),
                      fit_type = c("hal", "sl"),
                      sl_learners_density = NULL,
                      haldensify_args = list(
                        grid_type = "equal_range",
                        lambda_seq = exp(seq(-1, -13, length = 300))
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
  
  W <- W_g

  # make data objects from inputs
  data_in <- data.table::as.data.table(cbind(A, W))
  if (!is.matrix(W)) W <- as.matrix(W)
  data.table::setnames(data_in, c("A", colnames(W)))
  data.table::set(data_in, j = "ipc_weights", value = samp_weights)

  # need a data set with the exposure stochastically shifted DOWNWARDS A-delta
  data_in_downshifted <- data.table::copy(data_in)
  data.table::set(data_in_downshifted, j = "A", value = shift_additive(
    A = data_in$A, delta = -delta
  ))

  # need a data set with the exposure stochastically shifted UPWARDS A+delta
  data_in_upshifted <- data.table::copy(data_in)
  data.table::set(data_in_upshifted, j = "A", value = shift_additive(
    A = data_in$A, delta = delta
  ))

  # need a data set with the exposure stochastically shifted UPWARDS A+2delta
  data_in_upupshifted <- data.table::copy(data_in)
  data.table::set(data_in_upupshifted, j = "A", value = shift_additive(
    A = data_in$A, delta = 2 * delta
  ))

  # if fitting sl3 density make sl3 tasks from the data
  if (fit_type == "sl" & !is.null(sl_learners_density)) {
    # sl3 task for data with natural (unshifted) exposure
    sl_task <- sl3::sl3_Task$new(
      data = data_in,
      outcome = "A",
      covariates = colnames(W),
      weights = "ipc_weights"
    )

    # sl3 task for data with exposure shifted DOWNWARDS A-delta
    sl_task_downshifted <- sl3::sl3_Task$new(
      data = data_in_downshifted,
      outcome = "A",
      covariates = colnames(W),
      weights = "ipc_weights"
    )

    # sl3 task for data with exposure shifted UPWARDS A+delta
    sl_task_upshifted <- sl3::sl3_Task$new(
      data = data_in_upshifted,
      outcome = "A",
      covariates = colnames(W),
      weights = "ipc_weights"
    )

    # sl3 task for data with exposure shifted UPWARDS A+2delta
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
        wts = list(samp_weights),
        haldensify_args
      ),
      recursive = FALSE
    )
    fit_g_exp_dens_hal <- do.call(haldensify::haldensify, fit_args)
  } else if (fit_type == "sl" & !is.null(sl_learners_density)) {
    suppressMessages(
      fit_g_exp_dens_sl <- sl_learners_density$train(sl_task)
    )
  }

  # predict probabilities for the UNSHIFTED data (A = a)
  if (fit_type == "hal" & is.null(sl_learners_density)) {
    pred_g_exp_noshift <-
      stats::predict(
        object = fit_g_exp_dens_hal,
        new_A = as.numeric(data_in$A),
        new_W = as.matrix(data_in[, colnames(W), with = FALSE]),
        trim = FALSE
      )
  } else if (fit_type == "sl" & !is.null(sl_learners_density)) {
    suppressMessages(
      pred_g_exp_noshift <- fit_g_exp_dens_sl$predict()
    )
  }

  # predict probabilities for the DOWNSHIFTED data (A = a - delta)
  if (fit_type == "hal" & is.null(sl_learners_density)) {
    pred_g_exp_downshifted <-
      stats::predict(
        object = fit_g_exp_dens_hal,
        new_A = as.numeric(data_in_downshifted$A),
        new_W = as.matrix(data_in[, colnames(W), with = FALSE]),
        trim = FALSE
      )
  } else if (fit_type == "sl" & !is.null(sl_learners_density)) {
    suppressMessages(
      pred_g_exp_downshifted <- fit_g_exp_dens_sl$predict(sl_task_downshifted)
    )
  }

  # predict probabilities for the UPSHIFTED data (A = a + delta)
  if (fit_type == "hal" & is.null(sl_learners_density)) {
    pred_g_exp_upshifted <-
      stats::predict(
        object = fit_g_exp_dens_hal,
        new_A = as.numeric(data_in_upshifted$A),
        new_W = as.matrix(data_in[, colnames(W), with = FALSE]),
        trim = FALSE
      )
  } else if (fit_type == "sl" & !is.null(sl_learners_density)) {
    suppressMessages(
      pred_g_exp_upshifted <- fit_g_exp_dens_sl$predict(sl_task_upshifted)
    )
  }

  # predict probabilities for the UPSHIFTED data (A = a + 2*delta)
  if (fit_type == "hal" & is.null(sl_learners_density)) {
    pred_g_exp_upupshifted <-
      stats::predict(
        object = fit_g_exp_dens_hal,
        new_A = as.numeric(data_in_upupshifted$A),
        new_W = as.matrix(data_in[, colnames(W), with = FALSE]),
        trim = FALSE
      )
  } else if (fit_type == "sl" & !is.null(sl_learners_density)) {
    suppressMessages(
      pred_g_exp_upupshifted <- fit_g_exp_dens_sl$predict(sl_task_upupshifted)
    )
  }

  # create output data.tables
  out <- data.table::as.data.table(cbind(
    pred_g_exp_downshifted,
    pred_g_exp_noshift,
    pred_g_exp_upshifted,
    pred_g_exp_upupshifted
  ))
  data.table::setnames(out, c("downshift", "noshift", "upshift", "upupshift"))
  return(out)
}

###############################################################################

#' Estimate the Censoring Mechanism
#'
#' @details Compute the censoring mechanism for the observed data, in order to
#'  apply a joint intervention for removing censoring by re-weighting.
#'
#' @param C_cens A \code{numeric} vector of loss to follow-up indicators.
#' @param A A \code{numeric} vector of observed exposure values.
#' @param W A \code{numeric} matrix of observed baseline covariate values.
#' @param samp_weights A \code{numeric} vector of observation-level sampling
#'  weights, as produced by the internal procedure to estimate the two-phase
#'  sampling mechanism \code{\link{est_samp}}.
#' @param fit_type A \code{character} indicating whether to use GLMs or Super
#'  Learner to fit the censoring mechanism. If option "glm" is selected, the
#'  argument \code{glm_formula} must NOT be \code{NULL}, instead containing a
#'  model formula (as per \code{\link[stats]{glm}}) as a \code{character}. If
#'  the option "sl" is selected, the argument \code{sl_learners} must NOT be
#'  \code{NULL}; instead, an instantiated \pkg{sl3} \code{Lrnr_sl} object,
#'  specifying learners and a metalearner for the Super Learner fit, must be
#'  provided. Consult the documentation of \pkg{sl3} for details.
#' @param glm_formula A \code{character} giving a \code{\link[stats]{formula}}
#'  for fitting a (generalized) linear model via \code{\link[stats]{glm}}.
#' @param sl_learners Object containing a set of instantiated learners from the
#'  \pkg{sl3}, to be used in fitting an ensemble model.
#'
#' @importFrom stats glm as.formula predict
#' @importFrom data.table as.data.table setnames copy set
#' @importFrom stringr str_detect
#' @importFrom assertthat assert_that
#'
#' @return A \code{numeric} vector of the propensity score for censoring.
est_g_cens <- function(C_cens,
                       A,
                       W,
                       samp_weights = rep(1, length(C_cens)),
                       fit_type = c("sl", "glm"),
                       glm_formula = "C_cens ~ .",
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

  # generate the data objects for fitting the outcome regression
  data_in <- data.table::as.data.table(cbind(C_cens, A, W))
  if (!is.matrix(W)) W <- as.matrix(W)
  data.table::setnames(data_in, c(
    "C_cens", "A",
    paste0("W", seq_len(ncol(W)))
  ))
  names_W <- colnames(data_in)[stringr::str_detect(colnames(data_in), "W")]
  data.table::set(data_in, j = "ipc_weights", value = samp_weights)

  # fit a logistic regression and extract the predicted probabilities
  if (fit_type == "glm" & !is.null(glm_formula)) {
    # need to remove IPCW column from input data.table object for GLM fits
    data.table::set(data_in, j = "ipc_weights", value = NULL)

    # fit a logistic regression for the (scaled) outcome
    suppressWarnings(
      fit_g_cens <- stats::glm(
        formula = stats::as.formula(glm_formula),
        family = "binomial",
        data = data_in,
        weights = samp_weights
      )
    )

    # predict conditional censoring mechanism
    suppressWarnings(
      pred_g_cens <- unname(stats::predict(
        object = fit_g_cens,
        newdata = data_in,
        type = "response"
      ))
    )
  }

  # fit a binary Super Learner and get the predicted probabilities
  if (fit_type == "sl" & !is.null(sl_learners)) {
    # make sl3 task for original data
    task_cens <- sl3::sl3_Task$new(
      data = data_in,
      covariates = c("A", names_W),
      outcome = "C_cens",
      outcome_type = "binomial",
      weights = "ipc_weights"
    )

    # fit new Super Learner to the data and predict
    sl_fit_cens <- sl_learners$train(task_cens)
    pred_g_cens <- sl_fit_cens$predict()
  }

  # generate output
  return(pred_g_cens)
}

###############################################################################

#' Estimate the Outcome Mechanism
#'
#' @details Compute the outcome regression for the observed data, including
#'  with the shift imposed by the intervention. This returns the outcome
#'  regression for the observed data (at A) and under the counterfactual shift
#'  shift (at A + delta).
#'
#' @param Y A \code{numeric} vector of observed outcomes.
#' @param C_cens A \code{numeric} vector of loss to follow-up indicators.
#' @param A A \code{numeric} vector of observed exposure values.
#' @param W A \code{numeric} matrix of observed baseline covariate values.
#' @param delta A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the exposure \code{A}. This is passed to the internal
#'  \code{\link{shift_additive}} and is currently limited to additive shifts.
#' @param samp_weights A \code{numeric} vector of observation-level sampling
#'  weights, as produced by the internal procedure to estimate the two-phase
#'  sampling mechanism \code{\link{est_samp}}.
#' @param fit_type A \code{character} indicating whether to use GLMs or Super
#'  Learner to fit the outcome regression. If the option "glm" is selected, the
#'  argument \code{glm_formula} must NOT be \code{NULL}, instead containing a
#'  model formula (as per \code{\link[stats]{glm}}) as a \code{character}. If
#'  the option "sl" is selected, the argument \code{sl_learners} must NOT be
#'  \code{NULL}; instead, an instantiated \pkg{sl3} \code{Lrnr_sl} object,
#'  specifying learners and a metalearner for the Super Learner fit, must be
#'  provided. Consult the documentation of \pkg{sl3} for details.
#' @param glm_formula A \code{character} giving a \code{\link[stats]{formula}}
#'  for fitting a (generalized) linear model via \code{\link[stats]{glm}}.
#' @param sl_learners Object containing a set of instantiated learners from the
#'  \pkg{sl3}, to be used in fitting an ensemble model.
#' @param glm_family The family to be used for glm estimation of Q.
#'
#' @importFrom stats glm as.formula predict
#' @importFrom data.table as.data.table setnames copy set
#' @importFrom stringr str_detect
#' @importFrom assertthat assert_that
#'
#' @return A \code{data.table} with two columns, containing estimates of the
#'  outcome mechanism at the natural value of the exposure Q(A, W) and an
#'  upshift of the exposure Q(A + delta, W).
est_Q <- function(Y,
                  C_cens = rep(1, length(Y)),
                  A,
                  W,
                  delta = 0,
                  samp_weights = rep(1, length(Y)),
                  fit_type = c("sl", "glm"),
                  glm_formula = "Y ~ .",
                  glm_family = "binomial",
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
  # don't scale before prediction
  #y_star <- scale_to_unit(vals = Y)

  # generate the data objects for fitting the outcome regression
  #data_in <- data.table::as.data.table(cbind(y_star, C_cens, A, W))
  data_in <- data.table::as.data.table(cbind(Y, C_cens, A, W))
  if (!is.matrix(W)) W <- as.matrix(W)
  data.table::setnames(data_in, c(
    "Y", "C_cens", "A",
    paste0("W", seq_len(ncol(W)))
  ))
  names_W <- colnames(data_in)[stringr::str_detect(colnames(data_in), "W")]
  data.table::set(data_in, j = "ipc_weights", value = samp_weights)

  # counterfactual data with exposure shifted do(A + delta)
  # NOTE: also introduce joint intervention on the censoring node do(C = 1)
  data_in_shifted <- data.table::copy(data_in)
  data.table::set(data_in_shifted, j = "A", value = shift_additive(
    A = data_in$A, delta = delta
  ))
  data.table::set(data_in_shifted, j = "C_cens", value = 1)

  # fit a glm and extract the predicted probabilities
  if (fit_type == "glm" & !is.null(glm_formula)) {
    # need to remove IPCW column from input data.table object for GLM fits
    data.table::set(data_in, j = "ipc_weights", value = NULL)
    data.table::set(data_in_shifted, j = "ipc_weights", value = NULL)

    # fit a logistic regression for the (scaled) outcome
    suppressWarnings(
      fit_Qn <- stats::glm(
        formula = stats::as.formula(glm_formula),
        data = data_in,
        weights = samp_weights,
        subset = C_cens == 1,
        family = glm_family
      )
    )

    # joint intervention to remove censoring before prediction
    data.table::set(data_in, j = "C_cens", value = 1)

    # predict Qn for the un-shifted data do(A = a) & do(C = 1)
    suppressWarnings(
      pred_star_Qn <- unname(stats::predict(
        object = fit_Qn,
        newdata = data_in,
        type = "response"
      ))
    )

    suppressWarnings(
      # predict Qn for the shifted data (A = a + delta)
      pred_star_Qn_shifted <- unname(stats::predict(
        object = fit_Qn,
        newdata = data_in_shifted,
        type = "response"
      ))
    )
  }

  # fit a binary Super Learner and get the predicted probabilities
  if (fit_type == "sl" & !is.null(sl_learners)) {
    # make sl3 task for original data
    task_noshift <- sl3::sl3_Task$new(
      data = data_in[C_cens == 1, ],
      covariates = c("C_cens", "A", names_W),
      outcome = "Y",
      weights = "ipc_weights"
    )
    task_noshift_nocens <- sl3::sl3_Task$new(
      data = data_in[, C_cens := 1],
      covariates = c("C_cens", "A", names_W),
      outcome = "Y",
      weights = "ipc_weights"
    )

    # make sl3 task for data with the shifted exposure
    task_shifted <- sl3::make_sl3_Task(
      data = data_in_shifted,
      covariates = c("C_cens", "A", names_W),
      outcome = "Y",
      weights = "ipc_weights"
    )

    # fit new Super Learner to the natural (no shift) data and predict
    sl_fit_noshift <- sl_learners$train(task_noshift)
    sl_fit_noshift$print()
    pred_star_Qn <- sl_fit_noshift$predict(task_noshift_nocens)

    # predict with Super Learner from unshifted data on the shifted data
    pred_star_Qn_shifted <- sl_fit_noshift$predict(task_shifted)
  }
  
  #scale estimates
  pred_star_Qn <- (pred_star_Qn - min(Y)) / (max(Y) - min(Y))
  pred_star_Qn_shifted <- ( pred_star_Qn_shifted - min(Y)) / (max(Y) - min(Y))

  # NOTE: clean up predictions and generate output
  # avoid values that are exactly 0 or 1 in the scaled estimates
  pred_star_Qn <- bound_precision(vals = pred_star_Qn)
  pred_star_Qn_shifted <- bound_precision(vals = pred_star_Qn_shifted)

  # create output data frame and return result
  out <- data.table::as.data.table(cbind(pred_star_Qn, pred_star_Qn_shifted))
  data.table::setnames(out, c("noshift", "upshift"))
  return(out)
}

###############################################################################

#' Estimate Probability of Censoring by Two-Phase Sampling
#'
#' @details Compute estimates of the sampling probability for inclusion in the
#'  the second-phase via the two-phase sampling mechanism. These estimates are
#'  used for the creation of inverse probability weights.
#'
#' @param V A \code{numeric} vector, \code{matrix}, \code{data.frame} or
#'  similar object giving the observed values of the covariates known to
#'  potentially inform the sampling mechanism.
#' @param C_samp A \code{numeric} vector of observed values of the indicator
#'  for inclusion in the second-phase sample.
#' @param fit_type A \code{character} indicating whether to perform the fit
#'  using GLMs or a Super Learner ensemble model. If use of Super Learner is
#'  desired, then the argument \code{sl_learners} must be provided.
#' @param sl_learners An \pkg{sl3} \code{Lrnr_sl} object, a Super Learner
#'  ensemble or learner instantiated externally using \pkg{sl3}.
#'
#' @importFrom stats glm predict binomial
#' @importFrom data.table as.data.table setnames
#' @importFrom assertthat assert_that
#'
#' @return A \code{numeric} vector of the estimated sampling mechanism.
est_samp <- function(V,
                     C_samp,
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
  data_in <- data.table::as.data.table(cbind(C_samp, V))
  if (!is.matrix(V)) V <- as.matrix(V)
  names_V <- paste0("V", seq_len(ncol(V)))
  data.table::setnames(data_in, c("C_samp", names_V))

  # make sl3 tasks from the data if fitting Super Learners
  if (fit_type == "sl" & !is.null(sl_learners)) {
    sl_task <- sl3::sl3_Task$new(
      data = data_in,
      covariates = names_V,
      outcome = "C_samp",
      outcome_type = "binomial"
    )
  }

  # fit logistic regression or binomial SL to get class probabilities for IPCW
  if (fit_type == "glm" & is.null(sl_learners)) {
    samp_reg <- stats::glm(
      as.formula("C_samp ~ ."),
      family = stats::binomial(),
      data = data_in
    )
    samp_probs <- unname(stats::predict(
      object = samp_reg,
      newdata = data_in,
      type = "response"
    ))
  } else if (fit_type == "sl" & !is.null(sl_learners)) {
    sl_fit <- sl_learners$train(sl_task)
    samp_probs <- as.numeric(sl_fit$predict())
  }

  # the estimated sampling weights Pi_n for creating inverse weights
  return(samp_probs)
}

###############################################################################

#' Estimate Auxiliary Covariate of Full Data Efficient Influence Function
#'
#' @details Compute an estimate of the auxiliary covariate of the efficient
#'  influence function required to update initial estimates through logistic
#'  tilting models for targeted minimum loss estimation.
#'
#' @param gn_exp An estimate of the exposure density (a generalized propensity
#'  score) using the output provided by \code{\link{est_g_exp}}.
#' @param gps_bound \code{numeric} giving the lower limit of the generalized
#'  propensity score estimates to be tolerated (default = 0.05). Estimates less
#'  than this are truncated to this or 1/n. See \code{\link{bound_propensity}}
#'  for details.
#'
#' @importFrom data.table as.data.table setnames
#'
#' @return A \code{data.table} with two columns, containing estimates of the
#'  auxiliary covariate at the natural value of the exposure H(A, W) and at the
#'  shifted value of the exposure H(A + delta, W).
est_Hn <- function(gn_exp, gps_bound) {
  # set any g(a|w) = 0 values to a very small value above zero
  gn_exp$noshift <- bound_propensity(gn_exp$noshift, gps_bound)

  # compute the ratio of the propensity scores for Hn(A,W)
  ratio_g_noshift <- (gn_exp$downshift / gn_exp$noshift) +
    as.numeric(gn_exp$upshift == 0)

  # compute the ratio of the propensity scores for Hn(d(A,W),W)
  ratio_g_shift <- (gn_exp$noshift / gn_exp$upshift) *
    as.numeric(gn_exp$upshift != 0) + as.numeric(gn_exp$upupshift == 0)
  ratio_g_shift[is.na(ratio_g_shift)] <- 1

  # set up output table of auxiliary covariate under different shifts
  aux_covar <- data.table::as.data.table(cbind(ratio_g_noshift, ratio_g_shift))
  data.table::setnames(aux_covar, c("noshift", "shift"))
  return(aux_covar)
}
