#' Confidence Intervals for Shifted Treatment Parameters
#'
#' description
#'
#' @param object An object of class \code{shifttx}, as produced by invoking
#'  the function \code{tmle_shifttx}, for which a confidence interval is to be
#'  computed.
#' @param parm A \code{numeric} vector indicating indices of \code{object$est}
#'  for which to return confidence intervals.
#' @param level A \code{numeric} indicating the level of the confidence interval
#'  to be computed.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom stats qnorm
#'
#' @method confint shifttx
#'
#' @export
#
confint.shifttx <- function(object,
                            parm = seq_len(object$psi),
                            level = 0.95,
                            ...) {

    # first, let's get Z_(1 - alpha)
    norm_bounds <- c(-1, 1) * abs(stats::qnorm(p = (1 - level) / 2))

    # compute the EIF variance multiplier for the CI
    n_obs <- length(object$eif)
    sd_eif <- sqrt(object$var / n_obs)

    # compute the interval around the point estimate
    ci_psi <- norm_bounds * sd_eif + object$psi

    # set up output CI object
    ci_out <- c(ci_psi[1], object$psi, ci_psi[2])
    names(ci_out) <- c("lwr_CI", "est", "upr_CI")
    return(ci_out)
}

################################################################################

#' Print Method for shifttx Objects
#'
#' The \code{print} method for objects of class \code{shifttx}.
#'
#' @param x An object of class \code{shifttx}.
#' @param ... Other options (not currently used).
#'
#' @export
#'
#' @method print shifttx
#'

print.shifttx <- function(x, ...) {
  print(x[c("psi", "var")])
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

