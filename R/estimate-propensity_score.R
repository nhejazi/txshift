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
#' @importFrom data.table as.data.table setnames set copy
#' @importFrom condensier fit_density predict_probability speedglmR6
#' @importFrom sl3 sl3_Task
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

  # need a data set with the treatment stochastically shifted DOWNWARDS...
  data_in_downshifted <- data.table::copy(data_in)
  data.table::set(data_in_downshifted, j = "A", value = tx_shift(
    A = data_in$A, delta = delta, type = "additive", direc = "down"
  ))

  # need a data set with the treatment stochastically shifted UPWARDS...
  data_in_upshifted <- data.table::copy(data_in)
  data.table::set(data_in_upshifted, j = "A", value = tx_shift(
    A = data_in$A, delta = delta, type = "additive", direc = "up"
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

    # sl3 task for data with treatment shifted DOWNWARDS
    sl_task_downshifted <- sl3::sl3_Task$new(
      data = data_in_downshifted,
      outcome = "A",
      covariates = colnames(W),
      weights = "ipc_weights"
    )

    # sl3 task for data with treatment shifted UPWARDS
    sl_task_upshifted <- sl3::sl3_Task$new(
      data = data_in_upshifted,
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
  # create output data.tables: (1) A = a, (2) A = a - delta, (3) A = a + delta
  ##############################################################################
  out <- data.table::as.data.table(cbind(
    pred_g_A_downshifted,
    pred_g_A_noshift,
    pred_g_A_upshifted
  ))
  data.table::setnames(out, c("downshift", "noshift", "upshift"))
  return(out)
}
