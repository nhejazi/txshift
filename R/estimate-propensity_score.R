#' Estimate the Treatment Mechanism (Propensity Score)
#'
#' Compute the propensity score (treatment mechanism) for the observed data,
#' including the shift. This returns the propensity score for the observed data
#' (at A_i) and the shift of interest (at A_i - delta).
#'
#' @param A A \code{numeric} vector of observed treatment values.
#' @param W A \code{numeric} matrix of observed baseline covariate values.
#' @param delta A \code{numeric} ...
#' @param ipc_weights ...
#' @param fit_type A \code{character} specifying whether to use Super Learner to
#'  fit the density estimation.
#' @param std_args A \code{list} of arguments to be passed to \code{fit_density}
#'  from the \code{condensier} package when the argument \code{fit_type} is set
#'  to "glm" (not invoking Super Learner).
#'
#' @importFrom condensier fit_density predict_probability speedglmR6
#' @importFrom data.table as.data.table setnames set
#'
#' @keywords internal
#'
#' @export
#'
#' @author Nima Hejazi
#' @author David Benkeser
#
est_g <- function(A,
                  W,
                  delta = 0,
                  ipc_weights = rep(1, length(Y)),
                  fit_type = c("sl", "glm"),
                  sl_lrnrs = NULL,
                  sl_task = NULL,
                  std_args = list(nbins = 20,
                                  bin_method = "dhist",
                                  bin_estimator = condensier::speedglmR6$new(),
                                  parfit = FALSE)
                 ) {
  # make data object
  data_in <- data.table::as.data.table(cbind(A, W))
  if (!is.matrix(W)) W <- as.matrix(W)
  data.table::setnames(data_in, c("A", paste0("W", seq_len(ncol(W)))))
  data.table::set(data_in, j = "ipc_weights", value = ipc_weights)

  # if fitting sl3 density make sl3 task with data
  if (fit_type == "sl" & !is.null(sl_lrnrs) & is.null(sl_task)) {
    sl_task <- sl3::sl3_Task$new(
      data_in, outcome = "A",
      covariates = paste0("W", seq_len(ncol(W))),
      weights = "ipc_weights"
    )
  }

  # fit conditional density with condensier
  if (fit_type == "glm" & is.null(sl_lrnrs)) {
    fit_args <- unlist(list(X = c(paste0("W", seq_len(ncol(W)))),
                            Y = "A", std_args), recursive = FALSE)
    fit_args$input_data <- data_in
    fit_g_A <- do.call(condensier::fit_density, fit_args)
  } else if (fit_type == "sl" & !is.null(sl_lrnrs)) {
    sl_fit <- sl_lrnrs$train(task)
  }

  # predict probabilities for the un-shifted data (A = a)
  if (fit_type == "glm") {
    pred_g_A_noshift <-
      condensier::predict_probability(
        model_fit = fit_g_A,
        newdata = data_in
      )
  } else if (fit_type == "sl") {
    pred_g_A_noshift <- sl_fit$predict()
  }

  # predict probabilities for the DOWNSHIFTED data (A = a - delta)
  data_in_downshifted <- data.table::copy(data_in)
  data.table::set(data_in_downshifted, j = "A", value = tx_shift(
    A = data_in$A, delta = delta, type = "additive", direc = "down"
  ))
  ## create new SL task if required
  if (fit_type == "sl") {
    sl_task_downshifted <- sl3::sl3_Task$new(
      data_in_downshifted, outcome = "A",
      covariates = paste0("W", seq_len(ncol(W))),
      weights = "ipc_weights"
    )
  }

  # predict for the downshifted data
  if (fit_type == "glm") {
    pred_g_A_downshifted <-
      condensier::predict_probability(
        model_fit = fit_g_A,
        newdata = data_in_downshifted
      )
  } else if (fit_type == "sl") {
    pred_g_A_downshifted <- sl_fit$predict(sl_task_downshifted)
  }


  # predict probabilities for the UPSHIFTED data (A = a + delta)
  data_in_upshifted <- data.table::copy(data_in)
  data.table::set(data_in_upshifted, j = "A", value = tx_shift(
    A = data_in$A, delta = delta, type = "additive", direc = "up"
  ))
  ## create new SL task if required
  if (fit_type == "sl") {
    sl_task_upshifted <- sl3::sl3_Task$new(
      data_in_upshifted, outcome = "A",
      covariates = paste0("W", seq_len(ncol(W))),
      weights = "ipc_weights"
    )
  }

  # predict for the upshifted data
  if (fit_type == "glm") {
    pred_g_A_upshifted <-
      condensier::predict_probability(
        model_fit = fit_g_A,
        newdata = data_in_upshifted
      )
  } else if (fit_type == "sl") {
    pred_g_A_upshifted <- sl_fit$predict(sl_task_upshifted)
  }

  # create output matrix: scenarios A = a, A = a - delta
  out <- data.table::as.data.table(cbind(
    pred_g_A_downshifted,
    pred_g_A_noshift,
    pred_g_A_upshifted
  ))
  data.table::setnames(out, c("downshift", "noshift", "upshift"))
  return(out)
}
