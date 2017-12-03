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
#' @param fit_method A \code{character} indicating whether to use GLMs or Super
#'  Learner to fit the outcome regression. If the option "glm" is selected, the
#'  argument \code{glm_formula} must NOT be \code{NULL}, instead containing a
#'  model formula (in the style of \code{stats::glm}) as a \code{character}. If
#'  the option "sl" is selected, both of the arguments \code{sl_learners} and
#'  \code{sl_metalearner} must NOT be \code{NULL}, instead containing a set of
#'  learners and a metalearner for the Super Learner fit. Please consult the
#'  documentation of the \code{sl3} package for details on Super Learner fits.
#' @param glm_formula ...
#' @param sl_learners ...
#' @param sl_metalearner ...
#'
#' @importFrom stats glm as.formula predict
#' @importFrom sl3 make_sl3_Task make_learner Stack Lrnr_sl
#' @importFrom stringr str_detect
#' @importFrom data.table as.data.table setnames
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
                  ipc_weights = NULL,
                  fit_method = c("glm", "sl"),
                  glm_formula = "Y ~ .",
                  sl_learners = c(list("Lrnr_mean"),
                                  list("Lrnr_glm_fast")),
                  sl_metalearner = "nnls") {

    # scale the outcome for the logit transform
    y_star <- bound_scaling(Y = Y, scale = "zero_one")

    # make data object but using y_star rather than raw outcome
    data_O <- data.table::as.data.table(cbind(y_star, A, W))
    if (!is.matrix(W)) W <- as.matrix(W)
    data.table::setnames(data_O, c("Y", "A", paste0("W", seq_len(ncol(W)))))
    names_W <- colnames(data_O)[stringr::str_detect(colnames(data_O), "W")]

    # get the shifted treatment values
    a_shifted <- tx_shift(A = data_O$A, delta = delta,
                          type = "additive", direc = "up")

    # create a copy of the data for the shifted data set
    # and replace A with the shifted treatment (A = a + delta)
    data_O_shifted <- data_O
    data_O_shifted$A <- a_shifted

    if (fit_method == "glm" & !is.null(glm_formula)) {
        # obtain a logistic regression fit for the (scaled) outcome regression
        suppressWarnings(
          fit_Qn <- stats::glm(stats::as.formula(glm_formula),
                               data = data_O,
                               family = "binomial",
                               weights = ipc_weights)
        )

        # predict Qn for the un-shifted data (A = a)
        pred_star_Qn <- stats::predict(fit_Qn,
                                       newdata = data_O,
                                       type = "response")

        # predict Qn for the shifted data (A = a + delta)
        pred_star_Qn_shifted <- stats::predict(fit_Qn,
                                               newdata = data_O_shifted,
                                               type = "response")
    }

    if (fit_method == "sl" & !is.null(sl_learners) & !is.null(sl_metalearner)) {
        # add IPC weights to the data
        if (!is.null(ipc_weights)) {
          data_O$ipc_weights <- ipc_weights
          data_O_shifted$ipc_weights <- ipc_weights

          # make sl3 task for original data
          task_noshift <- sl3::make_sl3_Task(data = data_O,
                                             covariates = c("A", names_W),
                                             outcome = "Y",
                                             outcome_type = "quasibinomial",
                                             weights = "ipc_weights")

          # make sl3 task for data with the shifted treatment
          task_shifted <- sl3::make_sl3_Task(data = data_O_shifted,
                                             covariates = c("A", names_W),
                                             outcome = "Y",
                                             outcome_type = "quasibinomial",
                                             weights = "ipc_weights")
        } else {
          # make sl3 task for original data
          task_noshift <- sl3::make_sl3_Task(data = data_O,
                                             covariates = c("A", names_W),
                                             outcome = "Y",
                                             outcome_type = "quasibinomial")

          # make sl3 task for data with the shifted treatment
          task_shifted <- sl3::make_sl3_Task(data = data_O_shifted,
                                             covariates = c("A", names_W),
                                             outcome = "Y",
                                             outcome_type = "quasibinomial")
        }

        # create learners from arbitrary list and set up a stack
        # TODO: change to use sl3::make_learner_stack
        sl_lrnrs <- list()
        for (i in seq_along(sl_learners)) {
          sl_lrnrs[[i]] <- eval(parse(text = paste("sl3::Lrnr", sl_learners[i],
                                                   sep = "_")))
        }
        sl_lrnrs_ready <- lapply(sl_lrnrs, sl3::make_learner)
        stack <- sl3::make_learner(sl3::Stack, sl_lrnrs_ready)

        # extract meta-learner and create an sl3 Super Learner
        metalearner <- sl3::make_learner(eval(parse(text = paste("Lrnr",
                                                                 sl_metalearner,
                                                                 sep = "_"))))
        sl <- sl3::Lrnr_sl$new(learners = stack, metalearner = metalearner)

        # fit new Super Learner to the no-shift data and predict
        sl_fit_noshift <- sl$train(task_noshift)
        pred_star_Qn <- sl_fit_noshift$predict()

        # predict with Super Learner from unshifted data on shifted data
        pred_star_Qn_shifted <- sl_fit_noshift$predict(task_shifted)
    }

    # avoid values that are exactly 0 or 1 in the scaled Qn and Qn_shifted
    pred_star_Qn <- bound_precision(vals = as.numeric(pred_star_Qn))
    pred_star_Qn_shifted <- bound_precision(vals =
                                            as.numeric(pred_star_Qn_shifted))

    # create output data frame and return result
    out <- as.data.frame(cbind(pred_star_Qn, pred_star_Qn_shifted))
    colnames(out) <- c("noshift", "upshift")
    rownames(out) <- NULL
    return(out)
}

