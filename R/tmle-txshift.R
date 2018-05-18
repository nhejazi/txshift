#' Compute Targeted Minimum Loss-Based Estimate of Shifted Treatment Parameter
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
#' @param ipcw_fit_args A \code{list} of arguments, all but one of which are
#'  passed to \code{est_ipcw}. For details, please consult the documentation for
#'  \code{est_ipcw}. The first element of this (i.e., \code{fit_type}) is used
#'  to determine how this regression is fit: "glm" for generalized linear model,
#'  "sl" for a Super Learner, and "fit_spec" a user-specified input of the form
#'  produced by \code{est_ipcw}. NOTE THAT this first argument is not passed to
#'  \code{est_ipcw}.
#' @param g_fit_args A \code{list} of arguments, all but one of which are passed
#'  to \code{est_g}. For further details, please consult the documentation for
#'  \code{est_g}. The first element of this (i.e., \code{fit_type}) is used to
#'  determine how this regression is fit: "glm" for a generalized linear model
#'  for fitting densities via the \code{condensier} package, "sl" for \code{sl3}
#'  learners used to fit Super Learner to densities via \code{Lrnr_condensier},
#'  and "fit_spec" to user-specified input of the form produced by \code{est_g}.
#'  NOTE THAT this first argument is not passed to \code{est_g}.
#' @param Q_fit_args A \code{list} of arguments, all but one of which are passed
#'  to \code{est_Q}. For further details, please consult the documentation for
#'  \code{est_Q}. The first element of this (i.e., \code{fit_type}) is used to
#'  determine how this regression is fit: "glm" for a generalized linear model
#'  for the outcome regression, "sl" for \code{sl3} learners used to fit a Super
#'  Learner for the outcome regression, and "fit_spec" to user-specified input
#'  of the form produced by \code{est_Q}. NOTE THAT this first argument is not
#'  passed to \code{est_g}.
#' @param eif_reg_spec Whether a flexible function ought to be used in the
#'  computation of a targeting step for the censored data case. By default, the
#'  method used is a nonparametric regression based on the Highly Adaptive LASSO
#'  (from package \code{hal9001}) is used. If set to \code{FALSE}, then a simple
#'  linear regression model is assumed. In this step, the efficient influence
#'  function (EIF) is regressed against the covariates that contribute to the
#'  censoring mechanism (i.e., EIF ~ V | C = 1).
#' @param ipcw_efficiency Whether to invoke an augmentation of the IPCW-TMLE
#'  procedure that performs an iterative process to ensure efficiency of the
#'  resulting estimate. The default is \code{TRUE}; only set to \code{FALSE} if
#'  possible inefficiency of the IPCW-TMLE is not a concern.
#' @param ipcw_fit_spec User-specified version of the argument above for fitting
#'  the censoring mechanism (\code{ipcw_fit_args}). Consult the documentation
#'  for that argument for details on how to properly use this. In general, this
#'  should only be used by advanced users familiar with both the underlying
#'  theory and this software implementation of said theory.
#' @param gn_fit_spec User-specified version of the argument above for fitting
#'  the censoring mechanism (\code{g_fit_args}). Consult the documentation for
#'  that argument for details on how to properly use this. In general, this
#'  should only be used by advanced users familiar with both the underlying
#'  theory and this software implementation of said theory.
#' @param Qn_fit_spec User-specified version of the argument above for fitting
#'  the censoring mechanism (\code{Q_fit_args}). Consult the documentation for
#'  that argument for details on how to properly use this. In general, this
#'  should only be used by advanced users familiar with both the underlying
#'  theory and this software implementation of said theory.
#'
#' @importFrom condensier speedglmR6
#' @importFrom data.table as.data.table setnames
#' @importFrom stringr str_detect
#' @importFrom dplyr filter select "%>%"
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
                         eif_tol = 1 / length(Y),
                         max_iter = 1e4,
                         ipcw_fit_args = list(
                           fit_type = c("glm", "sl", "fit_spec"),
                           sl_lrnrs = NULL
                         ),
                         g_fit_args = list(
                           fit_type = c("glm", "sl", "fit_spec"),
                           nbins = 35,
                           bin_method = "dhist",
                           bin_estimator =
                             condensier::speedglmR6$new(),
                           parfit = FALSE,
                           sl_lrnrs_dens = NULL
                         ),
                         Q_fit_args = list(
                           fit_type = c("glm", "sl", "fit_spec"),
                           glm_formula = "Y ~ .",
                           sl_lrnrs = NULL
                         ),
                         eif_reg_spec = TRUE,
                         ipcw_efficiency = TRUE,
                         ipcw_fit_spec = NULL,
                         gn_fit_spec = NULL,
                         Qn_fit_spec = NULL) {
  ##############################################################################
  # TODO: check arguments and set up some objects for programmatic convenience
  ##############################################################################
  call <- match.call(expand.dots = TRUE)
  fluc_method <- match.arg(fluc_method)

  # dissociate fit type from other arguments to simplify passing to do.call
  ipcw_fit_type <- unlist(ipcw_fit_args[names(ipcw_fit_args) == "fit_type"],
    use.names = FALSE
  )
  g_fit_type <- unlist(g_fit_args[names(g_fit_args) == "fit_type"],
    use.names = FALSE
  )
  Q_fit_type <- unlist(Q_fit_args[names(Q_fit_args) == "fit_type"],
    use.names = FALSE
  )
  ipcw_fit_args <- ipcw_fit_args[names(ipcw_fit_args) != "fit_type"]
  g_fit_args <- g_fit_args[names(g_fit_args) != "fit_type"]
  Q_fit_args <- Q_fit_args[names(Q_fit_args) != "fit_type"]

  # coerce W to matrix and, if no names in W, assign them generically
  if (!is.matrix(W)) W <- as.matrix(W)
  W_names <- colnames(W)
  if (is.null(W_names)) {
    W_names <- paste0("W", seq_len(ncol(W)))
    colnames(W) <- W_names
  }

  ##############################################################################
  # perform sub-setting of data and implement IPC weighting if required
  ##############################################################################
  if (!all(C == 1) & !is.null(V)) {
    V_in <- data.table::as.data.table(mget(V))
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
    if (!is.null(ipcw_fit_spec) & ipcw_fit_type == "fit_spec") {
      ipcw_out <- ipcw_fit_spec
    } else {
      ipcw_out <- do.call(est_ipcw, ipcw_estim_args)
    }

    # extract IPC weights for censoring case and normalize weights
    cens_weights <- ipcw_out$ipc_weights
    cens_weights_norm <- cens_weights / sum(cens_weights)

    # remove column corresponding to indicator for censoring
    data_internal <- data.table::as.data.table(list(W, A = A, C = C, Y = Y)) %>%
      dplyr::filter(C == 1) %>%
      dplyr::select(-C) %>%
      data.table::as.data.table()
  } else {
    # if no censoring, we can just use IPC weights that are identically 1
    cens_weights <- C
    cens_weights_norm <- cens_weights / sum(cens_weights)
    data_internal <- data.table::as.data.table(list(W, A = A, Y = Y))
  }

  ##############################################################################
  # initial estimate of the treatment mechanism (propensity score)
  ##############################################################################
  if (!is.null(gn_fit_spec) & g_fit_type == "fit_spec") {
    gn_estim <- gn_fit_spec
  } else {
    gn_estim_in <- list(
      A = data_internal$A,
      W = data_internal[, W_names, with = FALSE],
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
    } else if (g_fit_type == "sl") {
      # if fitting SL, we can discard all the standard condensier-related args
      g_fit_args <- g_fit_args[stringr::str_detect(names(g_fit_args), "sl")]
      # reshapes list of args to make passing to do.call possible
      gn_estim_args <- unlist(list(gn_estim_in, g_fit_args), recursive = FALSE)
    }
    # pass the relevant args for computing the propensity score to do.call
    gn_estim <- do.call(est_g, gn_estim_args)
  }

  ##############################################################################
  # initial estimate of the outcome regression
  ##############################################################################
  if (!is.null(Qn_fit_spec) & Q_fit_type == "fit_spec") {
    Qn_estim <- Qn_fit_spec
  } else {
    Qn_estim_in <- list(
      Y = data_internal$Y,
      A = data_internal$A,
      W = data_internal[, W_names, with = FALSE],
      delta = delta,
      ipc_weights = cens_weights,
      fit_type = Q_fit_type
    )
    # reshape args to pass to the relevant function for the outcome regression
    Qn_estim_args <- unlist(list(Qn_estim_in, Q_fit_args), recursive = FALSE)
    # invoke function to estimate outcome regression via do.call
    Qn_estim <- do.call(est_Q, Qn_estim_args)
  }

  ##############################################################################
  # initial estimate of the auxiliary ("clever") covariate
  ##############################################################################
  Hn_estim <- est_Hn(gn = gn_estim)

  ##############################################################################
  # invoke efficient IPCW-TMLE, per Rose & van der Laan (2011), if necessary
  ##############################################################################
  n_steps <- 0  # define iteration outside to easily return in output object
  if (ipcw_efficiency & !all(C == 1) & !is.null(V)) {
    # Efficient implementation of the IPCW-TMLE
    eif_mean <- Inf
    conv_res <- replicate(3, rep(NA, max_iter))

    # normalize censoring mechanism weights (to be overwritten)
    cens_weights <- C / ipcw_out$pi_mech
    cens_weights_norm <- cens_weights / sum(cens_weights)

    # quantities to be updated in iterative procedure (to be overwritten)
    pi_mech_star <- ipcw_out$pi_mech
    Qn_estim_use <- Qn_estim

    # figure out columns of internal data structure used for censoring
    # TODO: FIX -- this is a horribly UGLY HACK that solves a naming problem...
    V_cols <- matrix(NA, nrow = ncol(V_in), ncol = ncol(data_internal))
    for (i in seq_along(V_in)) {
      V_cols[i, ] <- stringr::str_detect(
        colnames(V_in)[i],
        colnames(data_internal)
      )
    }
    V_mask <- as.logical(colSums(V_cols))
    colnames(V_in) <- colnames(data_internal[, which(V_mask), with = FALSE])

    # iterate procedure until convergence conditions are satisfied
    while (abs(eif_mean) > eif_tol & n_steps < max_iter) {
      # iterate counter
      n_steps <- n_steps + 1

      # update sub-model fluctuation, re-compute EIF, and update EIF
      ipcw_tmle_comp <- ipcw_tmle_proc(
        data_in = data_internal,
        C = C,
        V = V_in,
        ipc_mech = pi_mech_star,
        ipc_weights = cens_weights,
        ipc_weights_norm = cens_weights_norm,
        Qn_estim = Qn_estim_use,
        Hn_estim = Hn_estim,
        fluc_method = fluc_method,
        fit_type = ipcw_fit_type,
        eif_tol = eif_tol,
        sl_lrnrs = ipcw_fit_args$sl_lrnrs,
        eif_reg_spec = eif_reg_spec
      )

      # overwrite and update quantities to be used in the next iteration
      Qn_estim_use <- data.table::as.data.table(
        list(
          # NOTE: need to re-scale estimated outcomes values within bounds of Y
          bound_scaling(
            Y = Y,
            pred_vals = ipcw_tmle_comp$fluc_mod_out$Qn_noshift_star,
            scale_target =
              ipcw_tmle_comp$fluc_mod_out$Qn_noshift_star,
            scale_type = "bound_in_01"
          ),
          bound_scaling(
            Y = Y,
            pred_vals = ipcw_tmle_comp$fluc_mod_out$Qn_shift_star,
            scale_target =
              ipcw_tmle_comp$fluc_mod_out$Qn_shift_star,
            scale_type = "bound_in_01"
          )
        )
      )
      data.table::setnames(Qn_estim_use, names(Qn_estim))
      cens_weights <- ipcw_tmle_comp$ipc_weights
      cens_weights_norm <- ipcw_tmle_comp$ipc_weights_norm
      pi_mech_star <- ipcw_tmle_comp$pi_mech_star

      # compute updated mean of efficient influence function and save
      eif_mean <- mean(ipcw_tmle_comp$tmle_eif$eif - ipcw_tmle_comp$ipcw_eif)
      eif_var <- var(ipcw_tmle_comp$tmle_eif$eif - ipcw_tmle_comp$ipcw_eif) /
        length(C)
      conv_res[n_steps, ] <- c(ipcw_tmle_comp$tmle_eif$psi, eif_var, eif_mean)
    }
    conv_results <- data.table::as.data.table(conv_res)
    data.table::setnames(conv_results, c("psi", "var", "eif_mean"))
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
      fluc_mod_out = fitted_fluc_mod,
      Hn = Hn_estim,
      Y = data_internal$Y,
      Delta = C,
      ipc_weights = cens_weights,
      ipc_weights_norm = cens_weights_norm,
      tol_eif = eif_tol
    )
  }

  ##############################################################################
  # create output object
  ##############################################################################
  if (ipcw_efficiency & !all(C == 1) & !is.null(V)) {

    # replace variance in this object with the updated variance if iterative
    if (exists("eif_var")) {
      ipcw_tmle_comp$tmle_eif$var <- eif_var
    }

    # return only the useful convergence results
    conv_results_out <- conv_results[!is.na(rowSums(conv_results)), ]

    # create output object
    txshift_out <- unlist(
      list(
        call = call,
        ipcw_tmle_comp$tmle_eif,
        iter_res = list(conv_results_out),
        n_iter = n_steps
      ),
      recursive = FALSE
    )
  } else {
    # create output object
    txshift_out <- unlist(
      list(
        call = call,
        tmle_eif_out,
        n_iter = n_steps
      ),
      recursive = FALSE
    )
  }
  # S3-ify and return output object
  class(txshift_out) <- "txshift"
  return(txshift_out)
}
