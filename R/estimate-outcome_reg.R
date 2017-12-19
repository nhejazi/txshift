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
#' @param fit_type A \code{character} indicating whether to use GLMs or Super
#'  Learner to fit the outcome regression. If the option "glm" is selected, the
#'  argument \code{glm_formula} must NOT be \code{NULL}, instead containing a
#'  model formula (in the style of \code{stats::glm}) as a \code{character}. If
#'  the option "sl" is selected, the argument \code{sl_lrnrs} must NOT be
#'  \code{NULL}, instead being a pre-defined object of class \code{Lrnr_sl},
#'  specifying learners and a metalearner for the Super Learner fit. Consult the
#'  documentation of the \code{sl3} package for details on Super Learner fits.
#' @param glm_formula ...
#' @param sl_lrnrs ...
#'
#' @importFrom stats glm as.formula predict
#' @importFrom sl3 make_sl3_Task make_learner Stack Lrnr_sl
#' @importFrom stringr str_detect
#' @importFrom dplyr "%>%"
#' @importFrom data.table as.data.table setnames copy set
#'
#' @keywords internal
#'
#' @export
#'
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
  # scale the outcome for the logit transform
  y_star <- bound_scaling(Y = Y, scale = "zero_one")
  # generate the data objects for fitting the outcome regression
  data_in <- data.table::as.data.table(cbind(y_star, A, W))
  if (!is.matrix(W)) W <- as.matrix(W)
  data.table::setnames(data_in, c("Y", "A", paste0("W", seq_len(ncol(W)))))
  names_W <- colnames(data_in)[stringr::str_detect(colnames(data_in), "W")]
  data.table::set(data_in, j = "ipc_weights", value = ipc_weights)

  # need a data set with the treatment stochastically shifted UPWARDS...
  data_in_shifted <- data.table::copy(data_in)
  data.table::set(data_in_shifted, j = "A", value = tx_shift(
    A = data_in$A, delta = delta, type = "additive", direc = "up"
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
        stats::as.formula(glm_formula),
        data = data_in,
        family = "binomial",
        weights = ipc_weights
      )
    )

    # predict Qn for the un-shifted data (A = a)
    pred_star_Qn <- stats::predict(
      fit_Qn,
      newdata = data_in,
      type = "response"
    ) %>%
      as.numeric()

    # predict Qn for the shifted data (A = a + delta)
    pred_star_Qn_shifted <- stats::predict(
      fit_Qn,
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
    task_noshift <- sl3::make_sl3_Task(
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
