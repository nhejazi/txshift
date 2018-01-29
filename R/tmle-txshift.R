#' Compute Targeted Maximum Likelihood Estimate for Treatment Shift Parameter
#'
#' description THIS IS A USER-FACING WRAPPER FUNCTION
#'
#' @param W A \code{matrix} or \code{data.frame} corresponding to a set of
#'  baseline covariates.
#' @param A A \code{numeric} vector corresponding to a treatment variable. The
#'  parameter of interest is defined as a location shift of this quantity.
#' @param Y A \code{numeric} vector corresponding to an outcome variable.
#' @param C A \code{numeric} binary vector giving information on whether a given
#'  observation was subject to censoring. This is used to compute an IPCW-TMLE
#'  in cases where two-stage sampling is performed. The default assumes that no
#'  censoring was present (i.e., a two-stage design was NOT used). N.B., this is
#'  equivalent to the term %\Delta in the notation used in the original Rose and
#'  van der Laan manuscript that introduced/formulated IPCW-TML estimators.
#' @param V The covariates that are used in determining the sampling procedure
#'  that gives rise to censoring. The default is \code{NULL} and corresponds to
#'  scenarios in which there is no censoring (in which case all values in the
#'  preceding argument \code{C} must be uniquely 1. To specify this, pass in a
#'  NAMED \code{list} identifying variables amongst W, A, Y that are thought to
#'  have played a role in defining the sampling/censoring mechanism (C).
#' @param delta A \code{numeric} value indicating the shift in the treatment to
#'  be used in defining the target parameter. This is defined with respect to
#'  the scale of the treatment (A).
#' @param fluc_method The method to be used in submodel fluctuation step of
#'  the TMLE computation. The choices are "standard" and "weighted".
#' @param eif_tol A \code{numeric} giving the convergence criterion for the TML
#'  estimator. This is the the maximum mean of the efficient influence function
#'  (EIF) to be used in declaring convergence (theoretically, should be zero).
#' @param max_iter A \code{numeric} integer giving the maximum number of steps
#'  to be taken in iterating to a solution of the efficient influence function.
#' @param mod_args A nested \code{list} of \code{list}s that specifies arguments
#'  to be passed to the various internal functions for the estimation procedure.
#'  Each of sub-list corresponds to a single internal function. As such, for
#'  details on (1) \code{ipcw_fit}, see the documentation of \code{est_ipcw};
#'  (2) \code{g_fit}, see the documentation of \code{est_g}; (3) \code{Q_fit},
#'  see the documentation for \code{est_Q}.
#'
#' @importFrom condensier speedglmR6
#' @importFrom data.table as.data.table
#' @importFrom tibble as_tibble
#' @importFrom stringr str_detect
#' @importFrom dplyr filter select "%>%"
#' @importFrom rootSolve uniroot.all
#'
#' @return S3 object of class \code{txshift} containing the results of the
#'  procedure to compute a TML estimate of the treatment shift parameter.
#'
#' @export
#
tmle_txshift <- function(W,
                         A,
                         Y,
                         C = rep(1, length(Y)),
                         V = NULL,
                         delta = 0,
                         fluc_method = c("standard", "weighted"),
                         eif_tol = 1e-9,
                         max_iter = 1e4,
                         ipcw_fit_args = list(
                           fit_type = c("glm", "sl"),
                           glm_formula = "Delta ~ .",
                           sl_lrnrs = NULL
                         ),
                         g_fit_args = list(
                           fit_type = c("glm", "sl"),
                           nbins = 35,
                           bin_method = "dhist",
                           bin_estimator =
                             condensier::speedglmR6$new(),
                           parfit = FALSE,
                           sl_lrnrs_dens = NULL
                         ),
                         Q_fit_args = list(
                           fit_type = c("glm", "sl"),
                           glm_formula = "Y ~ .",
                           sl_lrnrs = NULL
                         )) {
  ##############################################################################
  # TODO: check arguments and set up some objects for programmatic convenience
  ##############################################################################
  call <- match.call(expand.dots = TRUE)
  fluc_method <- match.arg(fluc_method)

  # dissociate fit type from other arguments to simplify passing to do.call
  ipcw_fit_type <- unlist(ipcw_fit_args[names(ipcw_fit_args) == "fit_type"],
                          use.names = FALSE)
  g_fit_type <- unlist(g_fit_args[names(g_fit_args) == "fit_type"],
                       use.names = FALSE)
  Q_fit_type <- unlist(Q_fit_args[names(Q_fit_args) == "fit_type"],
                       use.names = FALSE)
  ipcw_fit_args <- ipcw_fit_args[names(ipcw_fit_args) != "fit_type"]
  g_fit_args <- g_fit_args[names(g_fit_args) != "fit_type"]
  Q_fit_args <- Q_fit_args[names(Q_fit_args) != "fit_type"]

  ##############################################################################
  # perform sub-setting of data and implement IPC weighting if required
  ##############################################################################
  if (any(unique(C) == 1) & !is.null(V)) {
    V_in <- data.table::as.data.table(V)
    ipcw_estim_in <- list(
      V = V_in, Delta = C,
      fit_type = ipcw_fit_type
    )
    # reshapes the list of args so that it can be passed to do.call
    ipcw_estim_args <- unlist(
      list(ipcw_estim_in, ipcw_fit_args),
      recursive = FALSE
    )
    # compute the IPC weights by passing all args to the relevant function
    ipcw_out <- do.call(est_ipcw, ipcw_estim_args)
    cens_weights <- ipcw_out$ipc_weights
    data_internal <- tibble::as_tibble(list(W = W, A = A, C = C, Y = Y)) %>%
      dplyr::filter(C == 1) %>%
      dplyr::select(-C) %>%
      data.table::as.data.table()
  } else {
    # if no censoring, we can just use IPC weights that are identically 1
    cens_weights <- C
    data_internal <- data.table::as.data.table(list(W = W, A = A, Y = Y))
  }

  ##############################################################################
  # initial estimate of the treatment mechanism (propensity score)
  ##############################################################################
  gn_estim_in <- list(
    A = data_internal$A,
    W = data_internal$W,
    delta = delta,
    ipc_weights = cens_weights,
    fit_type = g_fit_type
  )
  if (g_fit_type == "glm") {
    # since fitting a GLM, can safely remove all args related to SL
    g_fit_args <- g_fit_args[!stringr::str_detect(names(g_fit_args), "sl")]
    # reshape args to a list suitable to be passed to do.call
    gn_estim_args <- unlist(
      list(gn_estim_in, std_args = list(g_fit_args)),
      recursive = FALSE
    )
  } else {
    # if fitting SL, we can discard all the standard condensier-related args
    g_fit_args <- g_fit_args[stringr::str_detect(names(g_fit_args), "sl")]
    # reshapes list of args to make passing to do.call possible
    gn_estim_args <- unlist(list(gn_estim_in, g_fit_args), recursive = FALSE)
  }
  # pass the relevant args for computing the propensity score to do.call
  gn_estim <- do.call(est_g, gn_estim_args)

  ##############################################################################
  # initial estimate of the outcome regression
  ##############################################################################
  Qn_estim_in <- list(
    Y = data_internal$Y,
    A = data_internal$A,
    W = data_internal$W,
    delta = delta,
    ipc_weights = cens_weights,
    fit_type = Q_fit_type
  )
  # reshape args to pass to the relevant function for the outcome regression
  Qn_estim_args <- unlist(list(Qn_estim_in, Q_fit_args), recursive = FALSE)
  # invoke function to estimate outcome regression via do.call
  Qn_estim <- do.call(est_Q, Qn_estim_args)

  ##############################################################################
  # initial estimate of the auxiliary ("clever") covariate
  ##############################################################################
  Hn_estim <- est_Hn(gn = gn_estim)

  ##############################################################################
  # invoke efficient IPCW-TMLE, per Rose & van der Laan (2011), if necessary
  ##############################################################################
  if (any(unique(C) == 1) & !is.null(V)) {
    # Efficient implementation of the IPCW-TMLE
    n_steps <- 0
    eif_mean <- Inf
    conv_res <- rep(NA, max_iter)
    cens_weights_full <- C / ipcw_out$pi_mech

    # iterate procedure until convergence conditions are satisfied
    while (abs(eif_mean) > eif_tol & n_steps < max_iter) {
      # iterate counter
      n_steps <- n_steps + 1

      # update sub-model fluctuation, re-compute EIF, and update EIF
      ipcw_tmle_comp <- ipcw_tmle_proc(
        data_in = data_internal,
        C = C,
        V = V,
        ipcw_mech = ipcw_out$pi_mech,
        ipc_weights_all = cens_weights_full,
        Qn_estim = Qn_estim,
        Hn_estim = Hn_estim,
        fluc_method = fluc_method,
        fit_type = ipcw_fit_type,
        tol_eif = eif_tol,
        sl_lrnrs = ipcw_fit_args$sl_lrnrs
      )

      # compute updated mean of efficient influence function and save
      cens_weights_full <- ipcw_tmle_comp$ipc_weights
      eif_mean <- mean(ipcw_tmle_comp$tmle_eif$eif) - ipcw_tmle_comp$ipcw_eif
      conv_res[n_steps] <- eif_mean
    }
    ##############################################################################
    # standard TMLE of the shift parameter / inefficient IPCW-TMLE
    ##############################################################################
  } else {
    # fit logistic regression to fluctuate along the sub-model
    fitted_fluc_mod <- fit_fluc(
      Y = data_internal$Y,
      Qn_scaled = Qn_estim,
      Hn = Hn_estim,
      ipc_weights = cens_weights,
      method = fluc_method
    )
    # compute TML estimate and EIF for the treatment shift parameter
    tmle_eif_out <- tmle_eif(
      fluc_fit_out = fitted_fluc_mod,
      Hn = Hn_estim,
      Y = data_internal$Y,
      ipc_weights = cens_weights,
      tol_eif = eif_tol
    )
  }

  ##############################################################################
  # create output object
  ##############################################################################
  if (any(unique(C) == 1) & !is.null(V)) {
    # return only the useful convergence results
    conv_res_out <- conv_res[!is.na(conv_res)]
    txshift_out <- unlist(
      list(
        call = call,
        ipcw_tmle_comp$tmle_eif,
        eif_conv = list(conv_res_out)
      ),
      recursive = FALSE
    )
  } else {
    txshift_out <- unlist(list(call = call, tmle_eif_out), recursive = FALSE)
  }
  class(txshift_out) <- "txshift"
  return(txshift_out)
}

