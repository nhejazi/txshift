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
#' @param ...
#'
#' @importFrom condensier fit_density predict_probability
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
                  glm_formula = "A ~ .",
                  sl_lrnrs = NULL,
                  sl_task = NULL,
                  ...) {

  # make data object
  data_in <- data.table::as.data.table(cbind(A, W))
  if (!is.matrix(W)) W <- as.matrix(W)
  data.table::setnames(data_in, c("A", paste0("W", seq_len(ncol(W)))))

  # if fitting sl3 density make sl3 task with data
  if (fit_type == "sl" & !is.null(sl_lrnrs) & is.null(sl_task)) {
    data.table::set(data_in, j = "ipc_weights", ipc_weights)
    sl_task <- sl3::sl3_Task$new(
      data_in, outcome = "A",
      covariates = paste0("W", seq_len(ncol(W))),
      weights = "ipc_weights"
    )
  }

  # fit conditional density with condensier
  if (fit_type == "glm" & is.null(sl_lrnrs)) {
    fit_args <- unlist(list(X = c(paste0("W", seq_len(ncol(W)))),
                            Y = "A", input_data = data_in,
                            condensier_args), recursive = FALSE)
    fit_g_A <- do.call(condensier::fit_density, fit_args)
    fit_g_A <- condensier::fit_density(
      X = c(paste0("W", seq_len(ncol(W)))),
      Y = "A", input_data = data_in, ...
    )
  } else if (fit_type == "sl" & !is.null(lrnrs_sl)) {
    sl <- sl3::Lrnr_sl$new(
      learners = lrnrs_sl, metalearner =
        sl3::Lrnr_solnp_density$new()
    )
    sl_fit <- sl$train(task)
  }

  # predict probabilities for the un-shifted data (A = a)
  pred_g_A_noshift <- condensier::predict_probability(
    model_fit = fit_g_A,
    newdata = data_in
  )

  # predict probabilities for the DOWNSHIFTED data (A = a - delta)
  data_in_downshifted <- data.table::copy(data_in)
  data.table::set(data_in_downshifted, j = "A", value = tx_shift(
    A = data_in$A, delta = delta, type = "additive", direc = "down"
  ))
  pred_g_A_downshifted <-
    condensier::predict_probability(
      model_fit = fit_g_A,
      newdata = data_in_downshifted
    )

  # predict probabilities for the UPSHIFTED data (A = a + delta)
  data_in_upshifted <- data.table::copy(data_in)
  data.table::set(data_in_upshifted, j = "A", value = tx_shift(
    A = data_in$A, delta = delta, type = "additive", direc = "up"
  ))
  pred_g_A_upshifted <-
    condensier::predict_probability(
      model_fit = fit_g_A,
      newdata = data_in_upshifted
    )

  # create output matrix: scenarios A = a, A = a - delta
  out <- data.table::as.data.table(cbind(
    pred_g_A_downshifted,
    pred_g_A_noshift,
    pred_g_A_upshifted
  ))
  data.table::setnames(out, c("downshift", "noshift", "upshift"))
  return(out)
}
