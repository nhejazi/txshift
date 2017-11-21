#' Compute Targeted Minimum Loss-Based Estimates of Treatment Shift Parameters
#'
#' description THIS IS A USER-FACING WRAPPER FUNCTION
#'
#' @param W ...
#' @param A ...
#' @param Y ...
#' @param delta ...
#' @param parallel ...
#' @param Q_fit_method ...
#' @param fluc_method ...
#' @param eif_tol ...
#' @param ... Passed to internal function \code{est_g}, which then internally
#'  passes these arguments to \code{fit_density} from package \code{condensier}.
#'
#' @return S3 object of class \code{shifttx} containing the results of the
#'  procedure to compute a TML estimate of the treatment shift parameter.
#
tmle_shifttx <- function(W,
                         A,
                         Y,
                         delta,
                         parallel = NULL,
                         Q_fit_method = c("glm", "sl"),
                         fluc_method = c("standard", "weighted"),
                         eif_tol = 1e-7,
                         ...) {
    # check arguments
    # TODO

    # estimate the treatment mechanism (propensity score)
    gn_estim <- est_g(A = A, W = W, delta = delta, ...)

    # estimate the outcome regression
    Qn_estim <- est_Q(Y = Y, A = A, W = W, delta = delta,
                      fit_method = Q_fit_method)

    # estimate the auxiliary covariate
    Hn_estim <- est_Hn(gn = gn_estim)

    # fit the regression for submodel fluctuation
    fitted_fluc_mod <- fit_fluc(Y = Y,
                                Qn_scaled = Qn_estim,
                                Hn = Hn_estim,
                                method = fluc_method)

    # compute the TML estimate for the treatment shift parameter
    tml_estim <- tmle_eif(fluc_fit_out = fitted_fluc_mod, Hn = Hn_estim,
                          Y = Y, tol_eif = eif_tol)

    # create output object
    tmle_results <- list(tml_estim)
    class(tmle_results) <- "shifttx"
    return(tmle_results)
}

