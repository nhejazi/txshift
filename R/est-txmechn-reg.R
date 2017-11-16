#' Estimate the Treatment Mechanism
#'
#' Compute the propensity score (treatment mechanism) for the observed data,
#' including the shift. This returns the propensity score for the observed data
#' (at A_i) and the shift of interest (at A_i - delta).
#'
#' @param A A \code{numeric} vector of observed treatment values.
#' @param W A \code{numeric} matrix of observed baseline covariate values.
#' @param delta A \code{numeric} ...
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
                  ...) {
    # make data object
    data_O <- as.data.frame(cbind(A, W))
    if (!is.matrix(W)) W <- as.matrix(W)
    colnames(data_O) <- c("A", paste0("W", seq_len(ncol(W))))

    # fit conditional density with condensier
    fit_g_A <- fit_density(X = c(paste0("W", seq_len(ncol(W)))),
                           Y = "A", input_data = data_O, ...)

    # predict probabilities for the un-shifted data (A = a)
    pred_g_A_noshift <- predict_probability(model_fit = fit_g_A,
                                            newdata = data_O)

    # predict probabilities for the DOWNSHIFTED data (A = a - delta)
    data_O_downshifted <- data_O
    data_O_downshifted$A <- tx_shift(a = data_O_downshifted$A, delta = delta,
                                     type = "additive", direc = "down")
    pred_g_A_downshifted <- predict_probability(model_fit = fit_g_A,
                                                newdata = data_O_downshifted)

    # predict probabilities for the UPSHIFTED data (A = a + delta)
    data_O_upshifted <- data_O
    data_O_upshifted$A <- tx_shift(a = data_O_upshifted$A, delta = delta,
                                   type = "additive", direc = "up")
    pred_g_A_upshifted <- predict_probability(model_fit = fit_g_A,
                                              newdata = data_O_upshifted)

    # create output matrix: scenarios A = a, A = a - delta
    out <- as.data.frame(cbind(pred_g_A_downshifted, pred_g_A_noshift,
                               pred_g_A_upshifted))
    colnames(out) <- c("downshift", "noshift", "upshift")
    rownames(out) <- NULL
    return(out)
}

