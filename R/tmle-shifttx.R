#' Compute Targeted Maximum Likelihood Estimate for Treatment Shift Parameter
#'
#' description THIS IS A USER-FACING WRAPPER FUNCTION
#'
#' @param W ...
#' @param A ...
#' @param Y ...
#' @param C A \code{numeric} binary vector giving information on whether a given
#'  observation was subject to censoring. This is used to compute an IPCW-TMLE
#'  in cases where two-stage sampling is performed. The default assumes that no
#'  censoring was present (i.e., a two-stage design was NOT used). N.B., this is
#'  equivalent to the term %\Delta in the notation used in the original Rose and
#'  van der Laan manuscript that introduced/formulated IPCW-TML estimators.
#' @param delta ...
#' @param g_fit_args A \code{list} of arguments to be passed to the internal
#'  function \code{est_g}, which then passes these same arguments on to the
#'  function \code{fit_density} of the package \code{condensier}.
#' @param Q_fit_args A \code{list} of arugments to be passed to the function to
#'  estimate the outcome regression. Consult the documentation for the function
#'  \code{est_Q} for details.
#' @param fluc_method ...
#' @param eif_tol ...
#'
#' @importFrom dplyr filter "%>%"
#' @importFrom condensier speedglmR6
#' @importFrom data.table as.data.table
#'
#' @return S3 object of class \code{shifttx} containing the results of the
#'  procedure to compute a TML estimate of the treatment shift parameter.
#'
#' @export
#
tmle_shifttx <- function(W,
                         A,
                         Y,
                         delta,
                         C = rep(1, length(Y)),
                         ipcw_fit_args = list(
                           fit_type = "glm",
                           glm_formula = "Delta ~ .",
                           sl_lrnrs = NULL,
                           sl_task = NULL
                         ),
                         g_fit_args = list(
                           nbins = 20,
                           bin_method = "dhist",
                           bin_estimator =
                             condensier::speedglmR6$new(),
                           parfit = FALSE
                         ),
                         Q_fit_args = list(
                           fit_method = "glm",
                           glm_formula = "Y ~ .",
                           sl_lrnrs = NULL
                         ),
                         fluc_method = "standard",
                         eif_tol = 1e-7) {
  # TODO: check arguments

  # perform sub-setting of data and implement IPC weighting if required
  if (all(unique(C) != 1)) {
    ipcw_estim_in <- list(V = W, Delta = C)
    ipcw_estim_args <- unlist(list(ipcw_estim_in, ipcw_fit_args),
                              recursive = FALSE)
    cens_weights <- do.call(est_ipcw, ipcw_estim_args)
    O_nocensoring <- data.table::as.data.table(cbind(W, A, C, Y)) %>%
      dplyr::filter(C == 1)
  } else {
    cens_weights <- C
  }

  # estimate the treatment mechanism (propensity score)
  gn_estim_in <- list(A = A, W = W, delta = delta, ipc_weights = cens_weights)
  gn_estim_args <- unlist(list(gn_estim_in, g_fit_args), recursive = FALSE)
  gn_estim <- do.call(est_g, gn_estim_args)

  # estimate the outcome regression
  Qn_estim_in <- list(Y = Y, A = A, W = W, delta = delta)
  Qn_estim_args <- unlist(list(Qn_estim_in, Q_fit_args), recursive = FALSE)
  Qn_estim <- do.call(est_Q, Qn_estim_args)

  # estimate the auxiliary ("clever") covariate
  Hn_estim <- est_Hn(gn = gn_estim)

  # fit logistic regression to fluctuate along the sub-model
  fitted_fluc_mod <- fit_fluc(
    Y = Y,
    Qn_scaled = Qn_estim,
    Hn = Hn_estim,
    method = fluc_method
  )

  # compute Targeted Maximum Likelihood estimate for treatment shift parameter
  tmle_eif_out <- tmle_eif(
    Y = Y,
    Hn = Hn_estim,
    fluc_fit_out = fitted_fluc_mod,
    tol_eif = eif_tol
  )

  # create output object
  class(tmle_eif_out) <- "shifttx"
  return(tmle_eif_out)
}
