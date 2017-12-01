#' Fit Logistic Regression to Traverse the Fluctuation Submodel
#'
#' description
#'
#' @param Y ...
#' @param Qn_scaled ...
#' @param Hn ...
#' @param method ...
#'
#' @importFrom stats qlogis glm fitted predict
#'
#' @keywords internal
#
fit_fluc <- function(Y,
                     Qn_scaled,
                     Hn,
                     method = c("standard", "weighted")) {

    # scale the outcome for the logit transform
    y_star <- bound_scaling(Y = Y, scale = "zero_one")

    # extract Q and obtain logit transform
    logit_Qn <- stats::qlogis(Qn_scaled$noshift)

    # fit the fluctuation regression in one of two ways
    if (method == "standard") {
        # note that \epsilon_n will be the coefficient of the covariate Hn
        suppressWarnings(
          mod_fluc <- stats::glm(y_star ~ -1 + stats::offset(logit_Qn) +
                                   Hn$noshift,
                                 family = "binomial")
        )
    } else if (method == "weighted") {
        # note that \epsilon_n will be the intercept term here
        suppressWarnings(
          mod_fluc <- stats::glm(y_star ~ stats::offset(logit_Qn),
                                 weights = Hn$noshift,
                                 family = "binomial")
        )
    }

    # get fitted values from fluctuation model
    Qn_noshift_star_pred <- stats::fitted(mod_fluc)
    Qn_noshift_star <- bound_scaling(Y = Y,
                                     preds_scaled = Qn_noshift_star_pred,
                                     scale = "original")

    # need to logit transform Qn(d(A,W),W)
    Qn_shift_logit <- stats::qlogis(Qn_scaled$upshift)

    # get Qn_star for the SHIFTED data
    if (method == "standard") {
          Qn_shift_star_in <- as.data.frame(cbind(Qn_shift_logit, Hn$shift))
          colnames(Qn_shift_star_in) <- c("logit_Qn", "Hn$noshift")

          # predict from fluctuation model on Q(d(A,W),W) and Hn(d(A,W))
          Qn_shift_star_pred <- stats::predict(object = mod_fluc,
                                               newdata = Qn_shift_star_in,
                                               type = "response")

    } else if (method == "weighted") {
          Qn_shift_star_in <- as.data.frame(Qn_shift_logit)
          colnames(Qn_shift_star_in) <- "logit_Qn"

          # predict from fluctuation model on Q(d(A,W),W)
          Qn_shift_star_pred <- stats::predict(object = mod_fluc,
                                               newdata = Qn_shift_star_in,
                                               type = "response")
    }

    Qn_shift_star <- bound_scaling(Y = Y,
                                   preds_scaled = Qn_shift_star_pred,
                                   scale = "original")

    # return the fit model object
    out <- list(fluc_fit = mod_fluc,
                covar_method = method,
                Qn_shift_star = as.numeric(Qn_shift_star),
                Qn_noshift_star = as.numeric(Qn_noshift_star))
    return(out)
}

