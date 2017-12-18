#' Estimate Inverse Probability of Censoring Weights
#'
#' description
#'
#' @param V ...
#' @param Delta ...
#' @param fit_type A \code{character} indicating whether to perform the fit
#'  using GLMs or a Super Learner. If use of Super Learner is desired, then the
#'  argument \code{sl_lrnrs} must be provided.
#' @param sl_lrnrs An \code{R6} object of class \code{Lrnr_sl}, a Super Learner
#'  object created externally using the \code{sl3} package.
#' @param sl_task An \code{R6} object of class \code{sl3_Task} containing data
#'  to be used in fitting regressions via \code{sl3}. The default value of
#'  \code{NULL} should suffice in most cases as an appropriate \code{sl3_Task}
#'  is created automatically if the \code{sl_lrnrs} argument is NOT \code{NULL}.
#'
#' @importFrom stats glm predict binomial
#'
#' @return A \code{numeric} vector corresponding to the inverse probability of
#'  censoring weights required for computing an IPCW-TMLE. Formally, this is
#'  nothing more than %\frac{\Delta}{\Pi_n}, where the term %\Pi_n is simply the
#'  predicted probability of belonging to a censoring class as computed using
#'  standard logistic regression.
#'
#' @export
#
est_ipcw <- function(V,
                     Delta,
                     fit_type = c("sl", "glm"),
                     glm_formula = "Delta ~ .",
                     sl_lrnrs = NULL,
                     sl_task = NULL) {
  # coerce input V to matrix
  data_in <- data.table::as.data.table(cbind(Delta, V))
  if (!is.matrix(V)) V <- as.matrix(V)
  data.table::setnames(data_in, c("Delta", paste0("V", seq_len(ncol(V)))))

  # argument checks and setup for using SL
  if (fit_type == "sl" & !is.null(sl_lrnrs) & is.null(sl_task)) {
    sl_task <- sl3::sl3_Task$new(
      data_in, outcome = "Delta",
      covariates = paste0("V", seq_len(ncol(V))),
      outcome_type = "binomial"
    )
  }

  # fit logistic regression to get class probabilities for IPCW
  if (fit_type == "glm") {
    ipcw_reg <- stats::glm(as.formula(glm_formula),
                           family = stats::binomial(),
                           data = data_in)
    ipcw_probs <- stats::predict(
      object = ipcw_reg,
      newdata = data_in
    )
  } else if (fit_type == "sl" & !is.null(sl_lrnrs)) {
    sl_fit <- sl_lrnrs$train(sl_task)
    sl_fit_preds <- sl_fit$predict()
    ipcw_probs <- as.numeric(sl_fit_preds)
  } else {
    stop("Arguments for the model fitting process specified incorrectly.")
  }

  # compute the inverse weights as Delta/Pi_n and return this vector
  ipc_weights <- Delta / as.numeric(ipcw_probs)
  ipc_weights_out <- ipc_weights[ipc_weights != 0]
  return(ipc_weights_out)
}

