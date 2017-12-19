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
#' @param fit_type ...
#' @param fluc_method The method to be used in submodel fluctuation step of
#'  the TMLE computation. The choices are "standard" and "weighted".
#' @param eif_tol The convergence criterion for the TML estimator. This is the
#'  the maximum mean of the efficient influence function (EIF) to be accepted in
#'  calling convergence (theoretically, the value is zero).
#' @param args A nested \code{list} of \code{list}s that specifies arguments to
#'  be passed to the various internal functions for the estimation procedure.
#'  Each of sub-list corresponds to a single internal function. As such, for
#'  details on (1) \code{ipcw_fit}, see the documention of \code{est_ipcw}; (2)
#'  \code{g_fit}, see the documentation of \code{est_g}; (3) \code{Q_fit}, see
#'  the documentation for \code{est_Q}.
#'
#' @importFrom dplyr filter "%>%"
#' @importFrom condensier speedglmR6
#' @importFrom data.table as.data.table
#' @importFrom tibble as_tibble
#' @importFrom stringr str_detect
#'
#' @return S3 object of class \code{shifttx} containing the results of the
#'  procedure to compute a TML estimate of the treatment shift parameter.
#'
#' @export
#
tmle_shifttx <- function(W,
                         A,
                         Y,
                         C = rep(1, length(Y)),
                         delta = 0,
                         fit_type = c("glm", "sl"),
                         fluc_method = c("standard", "weighted"),
                         eif_tol = 1e-7,
                         args = list(
                           ipcw_fit = list(
                             glm_formula = "Delta ~ .",
                             sl_lrnrs = NULL
                           ),
                           g_fit = list(
                             nbins = 25,
                             bin_method = "dhist",
                             bin_estimator =
                               condensier::speedglmR6$new(),
                             parfit = FALSE,
                             sl_lrnrs_dens = NULL
                           ),
                           Q_fit = list(
                             glm_formula = "Y ~ .",
                             sl_lrnrs = NULL
                           )
                         )) {
  ##############################################################################
  # TODO: check arguments and set up some objects for programmatic convenience
  ##############################################################################
  call <- match.call(expand.dots = TRUE)
  fit_type <- match.arg(fit_type)
  fluc_method <- match.arg(fluc_method)

  # unpack the list of extra arguments for convenience
  ipcw_fit_args <- args$ipcw_fit
  g_fit_args <- args$g_fit
  Q_fit_args <- args$Q_fit

  ##############################################################################
  # perform sub-setting of data and implement IPC weighting if required
  ##############################################################################
  if (all(unique(C) != 1)) {
    ipcw_estim_in <- list(V = W, Delta = C, fit_type = fit_type)
    ipcw_estim_args <- unlist(
      list(ipcw_estim_in, ipcw_fit_args),
      recursive = FALSE
    )
    cens_weights <- do.call(est_ipcw, ipcw_estim_args)
    O_nocensoring <- tibble::as_tibble(list(W = W, A = A, C = C, Y = Y)) %>%
      dplyr::filter(C == 1)
  } else {
    cens_weights <- C
  }

  ##############################################################################
  # handle estimation under censoring with IPCW-TMLEs
  ##############################################################################
  if (all(unique(C) != 1)) {
    ############################################################################
    # estimate the treatment mechanism (propensity score)
    ############################################################################
    gn_estim_in <- list(
      A = O_nocensoring$A,
      W = O_nocensoring$W,
      delta = delta,
      ipc_weights = cens_weights,
      fit_type = fit_type
    )
    if (fit_type == "glm") {
      g_fit_args <- g_fit_args[!stringr::str_detect(names(g_fit_args), "sl")]
      gn_estim_args <- unlist(
        list(gn_estim_in, std_args = list(g_fit_args)),
        recursive = FALSE
      )
    } else {
      g_fit_args <- g_fit_args[stringr::str_detect(names(g_fit_args), "sl")]
      gn_estim_args <- unlist(list(gn_estim_in, g_fit_args), recursive = FALSE)
    }
    suppressMessages(
      gn_estim <- do.call(est_g, gn_estim_args)
    )

    ############################################################################
    # estimate the outcome regression
    ############################################################################
    Qn_estim_in <- list(
      Y = O_nocensoring$Y,
      A = O_nocensoring$A,
      W = O_nocensoring$W,
      delta = delta,
      fit_type = fit_type
    )
    Qn_estim_args <- unlist(list(Qn_estim_in, Q_fit_args), recursive = FALSE)
    Qn_estim <- do.call(est_Q, Qn_estim_args)

    ############################################################################
    # estimate the auxiliary ("clever") covariate
    ############################################################################
    Hn_estim <- est_Hn(gn = gn_estim)

    ############################################################################
    # fit logistic regression to fluctuate along the sub-model
    ############################################################################
    fitted_fluc_mod <- fit_fluc(
      Y = O_nocensoring$Y,
      Qn_scaled = Qn_estim,
      Hn = Hn_estim,
      method = fluc_method
    )

    ############################################################################
    # compute Targeted Maximum Likelihood estimate for treatment shift parameter
    ############################################################################
    tmle_eif_out <- tmle_eif(
      Y = O_nocensoring$Y,
      Hn = Hn_estim,
      fluc_fit_out = fitted_fluc_mod,
      tol_eif = eif_tol
    )
    ##############################################################################
    # when there is no censoring, just use a vanilla shift-TMLE
    ##############################################################################
  } else {
    ############################################################################
    # estimate the treatment mechanism (propensity score)
    ############################################################################
    gn_estim_in <- list(
      A = A,
      W = W,
      delta = delta,
      ipc_weights = cens_weights,
      fit_type = fit_type
    )
    if (fit_type == "glm") {
      g_fit_args <- g_fit_args[!stringr::str_detect(names(g_fit_args), "sl")]
      gn_estim_args <- unlist(
        list(gn_estim_in, std_args = list(g_fit_args)),
        recursive = FALSE
      )
    } else {
      g_fit_args <- g_fit_args[stringr::str_detect(names(g_fit_args), "sl")]
      gn_estim_args <- unlist(list(gn_estim_in, g_fit_args), recursive = FALSE)
    }
    suppressMessages(
      gn_estim <- do.call(est_g, gn_estim_args)
    )

    ############################################################################
    # estimate the outcome regression
    ############################################################################
    Qn_estim_in <- list(
      Y = Y,
      A = A,
      W = W,
      delta = delta,
      fit_type = fit_type
    )
    Qn_estim_args <- unlist(list(Qn_estim_in, Q_fit_args), recursive = FALSE)
    Qn_estim <- do.call(est_Q, Qn_estim_args)

    ############################################################################
    # estimate the auxiliary ("clever") covariate
    ############################################################################
    Hn_estim <- est_Hn(gn = gn_estim)

    ############################################################################
    # fit logistic regression to fluctuate along the sub-model
    ############################################################################
    fitted_fluc_mod <- fit_fluc(
      Y = Y,
      Qn_scaled = Qn_estim,
      Hn = Hn_estim,
      method = fluc_method
    )
    ############################################################################
    # compute Targeted Maximum Likelihood estimate for treatment shift parameter
    ############################################################################
    tmle_eif_out <- tmle_eif(
      Y = Y,
      Hn = Hn_estim,
      fluc_fit_out = fitted_fluc_mod,
      tol_eif = eif_tol
    )
  }
  ##############################################################################
  # create output object
  ##############################################################################
  shifttx_out <- unlist(list(call = call, tmle_eif_out), recursive = FALSE)
  class(shifttx_out) <- "shifttx"
  return(shifttx_out)
}
