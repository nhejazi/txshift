#' Estimate Inverse Probability of Censoring Weights
#'
#' description
#'
#' @param V ...
#' @param Delta ...
#' @param fit_type A \code{character} indicating whether to perform the fit
#'  using GLMs or a Super Learner. If use of Super Learner is desired, then the
#'  argument \code{sl_lrnrs} must be provided.
#' @param glm_formula ...
#' @param sl_lrnrs An \code{R6} object of class \code{Lrnr_sl}, a Super Learner
#'  object created externally using the \code{sl3} package.
#'
#' @importFrom stats glm predict binomial
#' @importFrom data.table as.data.table setnames
#' @importFrom dplyr "%>%"
#' @importFrom sl3 sl3_Task
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
                     sl_lrnrs = NULL) {
  ##############################################################################
  # make data objects from inputs
  ##############################################################################
  data_in <- data.table::as.data.table(cbind(Delta, V))
  if (!is.matrix(V)) V <- as.matrix(V)
  names_V <- paste0("V", seq_len(ncol(V)))
  data.table::setnames(data_in, c("Delta", names_V))

  ##############################################################################
  # make sl3 tasks from the data if fitting Super Learners
  ##############################################################################
  if (fit_type == "sl" & !is.null(sl_lrnrs)) {
    sl_task <- sl3::sl3_Task$new(
      data_in, outcome = "Delta",
      covariates = names_V,
      outcome_type = "binomial"
    )
  }

  ##############################################################################
  # fit logistic regression or binomial SL to get class probabilities for IPCW
  ##############################################################################
  if (fit_type == "glm" & is.null(sl_lrnrs)) {
    ipcw_reg <- stats::glm(
      as.formula(glm_formula),
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
  return(ipc_weights_out)
}

