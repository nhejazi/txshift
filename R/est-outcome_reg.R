#' Estimate the Outcome Regression
#'
#' Compute the outcome regression for the observed data, including with the
#' shift imposed by the intervention. This returns the propensity score for the
#' observed data (at A_i) and the shift of interest (at A_i - delta).
#'
#' @param Y A \code{numeric} vector of observed outcomes.
#' @param A A \code{numeric} vector of observed treatment values.
#' @param W A \code{numeric} matrix of observed baseline covariate values.
#' @param delta ...
#' @param glm_form ...
#' @param sl_learners ...
#' @param sl_meta ...
#'
#' @importFrom stats glm
#'
#' @keywords internal
#'
#' @export
#'
#' @author Nima Hejazi
#' @author David Benkeser
#
est_Q <- function(Y,
                  A,
                  W,
                  delta = 0,
                  glm_form = "Y ~ .",
                  sl_lrnrs = NULL,
                  sl_meta = NULL) {

    # scale the outcome for the logit transform
    y_star <- bound_scaling(Y = Y, scale = "zero_one")

    # make data object but using y_star rather than raw outcome
    data_O <- as.data.frame(cbind(y_star, A, W))
    if (!is.matrix(W)) W <- as.matrix(W)
    colnames(data_O) <- c("Y", "A", paste0("W", seq_len(ncol(W))))

    # get the shifted treatment values
    a_shifted <- tx_shift(a = data_O$A, delta = delta,
                          type = "additive", direc = "up")

    # create a copy of the data for the shifted data set
    # and replace A with the shifted treatment (A = a + delta)
    data_O_shifted <- data_O
    data_O_shifted$A <- a_shifted

    if (!is.null(glm_form)) {
        # obtain a logistic regression fit for the (scaled) outcome regression
        suppressWarnings(
          fit_Qn <- glm(as.formula(glm_form),
                        data = data_O,
                        family = "binomial")
        )

        # predict Qn for the un-shifted data (A = a)
        pred_star_Qn <- predict(fit_Qn,
                                newdata = data_O,
                                type = "response")

        # predict Qn for the shifted data (A = a + delta)
        pred_star_Qn_shifted <- predict(fit_Qn,
                                        newdata = data_O_shifted,
                                        type = "response")
    }

    if (!is.null(sl_lrnrs) & !is.null(sl_meta)) {
        # make sl3 task for original data
        task_noshift <- make_sl3_Task(data = data_O,
                                      covariates = c("A", paste0("W", seq_len(n_w))),
                                      outcome = "Y", outcome_type = "continuous")

        # make sl3 task for data with the shifted treatment
        task_shifted <- make_sl3_Task(data = data_O_shifted,
                                      covariates = c("A", paste0("W", seq_len(n_w))),
                                      outcome = "Y", outcome_type = "continuous")
        # fit SL
    }

    # avoid values that are exactly 0 or 1 in the scaled Qn and Qn_shifted
    pred_star_Qn <- bound_precision(vals = pred_star_Qn)
    pred_star_Qn_shifted <- bound_precision(vals = pred_star_Qn_shifted)

    # create output matrix: scenarios A = a, A = a - delta
    out <- as.data.frame(cbind(pred_star_Qn, pred_star_Qn_shifted))
    colnames(out) <- c("noshift", "upshift")
    rownames(out) <- NULL
    return(out)
}

