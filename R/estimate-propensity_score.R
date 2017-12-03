#' Estimate the Treatment Mechanism (Propensity Score)
#'
#' Compute the propensity score (treatment mechanism) for the observed data,
#' including the shift. This returns the propensity score for the observed data
#' (at A_i) and the shift of interest (at A_i - delta).
#'
#' @param A A \code{numeric} vector of observed treatment values.
#' @param W A \code{numeric} matrix of observed baseline covariate values.
#' @param delta A \code{numeric} ...
#' @param fit_type A \code{character} specifying whether to use Super Learner to
#'  fit the density estimation.
#' @param ...
#'
#' @importFrom condensier fit_density predict_probability
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
                  ipc_weights = NULL,
                  fit_type = "standard",
                  lrnrs_sl = NULL,
                  ...) {

    # make data object
    data_O <- data.table::as.data.table(cbind(A, W))
    if (!is.matrix(W)) W <- as.matrix(W)
    data.table::setnames(data_O, c("A", paste0("W", seq_len(ncol(W)))))

    # if fitting sl3 density make sl3 task with data
    if (fit_type == "sl") {
        if(!is.null(ipc_weights)) {
          data_O$ipc_weights <- ipc_weights
          task <- sl3::sl3_Task$new(data_O, outcome = "A",
                                    covariates = paste0("W", seq_len(ncol(W))),
                                    weights = "ipc_weights")
        } else {
          task <- sl3::sl3_Task$new(data_O, outcome = "A",
                                    covariates = paste0("W", seq_len(ncol(W))))
        }
    }

    # fit conditional density with condensier
    if (fit_type == "standard" & is.null(lrnrs_sl)) {
        fit_g_A <- condensier::fit_density(X = c(paste0("W", seq_len(ncol(W)))),
                                           Y = "A", input_data = data_O, ...)
    } else if (fit_type == "sl" & !is.null(lrnrs_sl)) {
        sl <- sl3::Lrnr_sl$new(learners = lrnrs_sl, metalearner =
                               sl3::Lrnr_solnp_density$new())
        sl_fit <- sl$train(task)
    }

    # predict probabilities for the un-shifted data (A = a)
    pred_g_A_noshift <- condensier::predict_probability(model_fit = fit_g_A,
                                                        newdata = data_O)

    # predict probabilities for the DOWNSHIFTED data (A = a - delta)
    data_O_downshifted <- data_O
    data_O_downshifted$A <- tx_shift(A = data_O_downshifted$A, delta = delta,
                                     type = "additive", direc = "down")
    pred_g_A_downshifted <-
        condensier::predict_probability(model_fit = fit_g_A,
                                        newdata = data_O_downshifted)

    # predict probabilities for the UPSHIFTED data (A = a + delta)
    data_O_upshifted <- data_O
    data_O_upshifted$A <- tx_shift(A = data_O_upshifted$A, delta = delta,
                                   type = "additive", direc = "up")
    pred_g_A_upshifted <-
        condensier::predict_probability(model_fit = fit_g_A,
                                        newdata = data_O_upshifted)

    # create output matrix: scenarios A = a, A = a - delta
    out <- as.data.frame(cbind(pred_g_A_downshifted, pred_g_A_noshift,
                               pred_g_A_upshifted))
    colnames(out) <- c("downshift", "noshift", "upshift")
    rownames(out) <- NULL
    return(out)
}

