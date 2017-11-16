#' Confidence Intervals for Shifted Treatment Parameters
#'
#' description
#'
#' @param tmle_shifttx_out ...
#' @param level ...
#'
#' @importFrom stats qnorm
#'
#' @keywords internal
#
confint.shifttx <- function(tmle_shifttx_out,
                            level = 0.95) {
    # first, let's get Z_(1 - alpha)
    norm_bounds <- c(-1, 1) * abs(stats::qnorm(p = (1 - level) / 2))

    # compute the EIF variance multiplier for the CI
    n_obs <- length(tmle_shifttx_out$eif)
    sd_eif <- sqrt(tmle_shifttx_out$var / n_obs)

    # compute the interval around the point estimate
    ci_psi <- norm_bounds * sd_eif + tmle_shifttx_out$psi

    # set up output CI object
    ci_out <- c(ci_psi[1], tmle_shifttx_out$psi, ci_psi[2])
    names(ci_out) <- c("lwrCI_psi", "est_psi", "uprCI_psi")
    return(ci_out)
}

################################################################################

#' Bound Precision
#'
#' description
#'
#' @param vals ...
#'
#' @keywords internal
#
bound_precision <- function(vals) {
    if (max(vals) > 1 | min(vals) < 0) {
        stop("Scaled values are not in the interval [0, 1].")
    }
    vals[vals == 0] <- .Machine$double.neg.eps
    vals[vals == 1] <- 1 - .Machine$double.neg.eps
    return(vals)
}

################################################################################

#' Induce Scaling
#'
#' description
#'
#' @param Y ...
#' @param preds_scaled ...
#' @param scale ...
#'
#' @keywords internal
#
bound_scaling <- function(Y,
                          preds_scaled = NULL,
                          scale = c("zero_one", "original")) {
    y_min <- min(Y)
    y_max <- max(Y)

    if (scale == "zero_one") {
        y_star <- (Y - y_min) / (y_max - y_min)
        return(y_star)
    } else if (scale == "original" & !is.null(preds_scaled)) {
        preds_original <- (y_max - y_min) * preds_scaled + y_min
        return(preds_original)
    }
}

