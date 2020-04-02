#' Estimate Counterfactual Mean Under Stochastic Shift in Exposure
#'
#' @details TODO
#'
#' @param W A \code{matrix}, \code{data.frame}, or similar containing a set of
#'  baseline covariates.
#' @param A A \code{numeric} vector corresponding to a treatment variable. The
#'  parameter of interest is defined as a location shift of this quantity.
#' @param Y A \code{numeric} vector of the observed outcomes.
#' @param C A \code{numeric} indicator for whether a given observation was
#'  subject to censoring, used to compute an IPC-weighted estimator in cases
#'  where two-stage sampling is performed. The default assumes no censoring.
#' @param V The covariates that are used in determining the sampling procedure
#'  that gives rise to censoring. The default is \code{NULL} and corresponds to
#'  scenarios in which there is no censoring (in which case all values in the
#'  preceding argument \code{C} must be uniquely 1). To specify this, pass in a
#'  \code{character} vector identifying variables amongst W, A, Y thought to
#'  have played a role in defining the sampling/censoring mechanism (C). This
#'  argument also accepts a \code{data.table} (or similar) object composed of
#'  combinations of variables W, A, Y; use of this option is NOT recommended.
#' @param delta A \code{numeric} value indicating the shift in the treatment to
#'  be used in defining the target parameter. This is defined with respect to
#'  the scale of the treatment (A).
#' @param estimator The type of estimator to be fit, either \code{"tmle"} for
#'  targeted maximum likelihood or \code{"onestep"} for a one-step estimator.
#' @param fluctuation The method to be used in the submodel fluctuation step
#'  (targeting step) to compute the TML estimator. The choices are "standard"
#'  and "weighted" for where to place the auxiliary covariate in the logistic
#'  tilting regression.
#' @param max_iter A \code{numeric} integer giving the maximum number of steps
#'  to be taken in iterating to a solution of the efficient influence function.
#' @param ipcw_fit_args A \code{list} of arguments, all but one of which are
#'  passed to \code{\link{est_ipcw}}. For details, consult the documentation of
#'  \code{\link{est_ipcw}}. The first element (i.e., \code{fit_type}) is used
#'  to determine how this regression is fit: generalized linear model ("glm")
#'  or Super Learner ("sl"), and "fit_spec" a user-specified input of the form
#'  produced by \code{\link{est_ipcw}}. NOTE THAT this first argument is not
#'  passed to \code{\link{est_ipcw}}.
#' @param g_fit_args A \code{list} of arguments, all but one of which are
#'  passed to \code{\link{est_g}}. For details, consult the documentation of
#'  \code{\link{est_g}}. The first element (i.e., \code{fit_type}) is used to
#'  determine how this regression is fit: \code{"hal"} to estimate conditional
#'  densities via the highly adaptive lasso (via \pkg{haldensify}), \code{"sl"}
#'  for \pkg{sl3} learners used to fit Super Learner to densities via
#'  \code{\link[sl3]{Lrnr_haldensify}} or similar, and \code{"fit_spec"} for
#'  user-specified input of the form produced by \code{\link{est_g}}. NOTE that
#'  this first argument is not passed to \code{\link{est_g}}.
#' @param Q_fit_args A \code{list} of arguments, all but one of which are
#'  passed to \code{\link{est_Q}}. For details, consult the documentation for
#'  \code{\link{est_Q}}. The first element (i.e., \code{fit_type}) is used to
#'  determine how this regression is fit: \code{"glm"} for a generalized linear
#'  model for the outcome regression, \code{"sl"} for \pkg{sl3} learners used
#'  to fit a Super Learner for the outcome regression, and \code{"fit_spec"}
#'  for user-specified input of the form produced by \code{\link{est_Q}}. NOTE
#'  that this first argument is not passed to \code{\link{est_g}}.
#' @param eif_reg_type Whether a flexible nonparametric function ought to be
#'  used in the dimension-reduced nuisance regression of the targeting step for
#'  the censored data case. By default, the method used is a nonparametric
#'  regression based on the Highly Adaptive Lasso (from \pkg{hal9001}). Set
#'  this to \code{"glm"} to instead use a simple linear regression model. In
#'  this step, the efficient influence function (EIF) is regressed against
#'  covariates contributing to the censoring mechanism (i.e., EIF ~ V | C = 1).
#' @param ipcw_efficiency Whether to invoke an augmentation of the IPCW-TMLE
#'  procedure that performs an iterative process to ensure efficiency of the
#'  resulting estimate. The default is \code{TRUE}; only set to \code{FALSE} if
#'  possible inefficiency of the IPCW-TMLE is not a concern.
#' @param ipcw_fit_spec User-specified version of the argument above for
#'  fitting the censoring mechanism (\code{ipcw_fit_args}). Consult the
#'  documentation for that argument for details on how to properly use this. In
#'  general, this should only be used by advanced users familiar with both the
#'  underlying theory and this software implementation.
#' @param gn_fit_spec User-specified version of the argument above for fitting
#'  the treatment mechanism (\code{g_fit_args}). Consult the documentation for
#'  that argument for details on how to properly use this. In general, this
#'  should only be used by advanced users familiar with both the underlying
#'  theory and this software implementation.
#' @param Qn_fit_spec User-specified version of the argument above for fitting
#'  the outcome mechanism (\code{Q_fit_args}). Consult the documentation for
#'  that argument for details on how to properly use this. In general, this
#'  should only be used by advanced users familiar with both the underlying
#'  theory and this software implementation.
#'
#' @importFrom data.table data.table as.data.table setnames ":="
#' @importFrom stringr str_detect
#' @importFrom Rdpack reprompt
#'
#' @return S3 object of class \code{txshift} containing the results of the
#'  procedure to compute a TML or one-step estimate of the counterfactual mean
#'  under a modified treatment policy that shifts a continuous-valued exposure
#'  by a scalar amount \code{delta}. These estimates can be augmented to be
#'  consistent and efficient when two-phase sampling is performed.
#'
#' @examples
#' set.seed(429153)
#' n_obs <- 100
#' W <- replicate(2, rbinom(n_obs, 1, 0.5))
#' A <- rnorm(n_obs, mean = 2 * W, sd = 1)
#' Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
#' C <- rbinom(n_obs, 1, plogis(W + Y)) # two-phase sampling
#'
#' # construct a TML estimate (set estimator = "onestep" for the one-step)
#' tmle <- txshift(
#'   W = W, A = A, Y = Y, delta = 0.5,
#'   estimator = "tmle",
#'   g_fit_args = list(
#'     fit_type = "hal", n_bins = 5,
#'     grid_type = "equal_mass",
#'     lambda_seq = exp(-1:-9)
#'   ),
#'   Q_fit_args = list(
#'     fit_type = "glm",
#'     glm_formula = "Y ~ ."
#'   )
#' )
#'
#' # construct a TML estimate under two-phase sampling
#' ipcwtmle <- txshift(
#'   W = W, A = A, Y = Y, delta = 0.5,
#'   C = C, V = c("W", "Y"),
#'   estimator = "tmle", max_iter = 5,
#'   ipcw_fit_args = list(fit_type = "glm"),
#'   g_fit_args = list(
#'     fit_type = "hal", n_bins = 5,
#'     grid_type = "equal_mass",
#'     lambda_seq = exp(-1:-9)
#'   ),
#'   Q_fit_args = list(
#'     fit_type = "glm",
#'     glm_formula = "Y ~ ."
#'   ),
#'   eif_reg_type = "glm"
#' )
#' @export
txshift <- function(W,
                    A,
                    Y,
                    C = rep(1, length(Y)),
                    V = NULL,
                    delta = 0,
                    estimator = c("tmle", "onestep"),
                    fluctuation = c("standard", "weighted"),
                    max_iter = 1e3,
                    ipcw_fit_args = list(
                      fit_type = c("glm", "sl", "fit_spec"),
                      sl_learners = NULL
                    ),
                    g_fit_args = list(
                      fit_type = c("hal", "sl", "fit_spec"),
                      n_bins = c(10, 25),
                      grid_type = c("equal_range", "equal_mass"),
                      lambda_seq = exp(seq(-1, -13, length = 300)),
                      use_future = FALSE,
                      sl_learners_density = NULL
                    ),
                    Q_fit_args = list(
                      fit_type = c("glm", "sl", "fit_spec"),
                      glm_formula = "Y ~ .",
                      sl_learners = NULL
                    ),
                    eif_reg_type = c("hal", "glm"),
                    ipcw_efficiency = TRUE,
                    ipcw_fit_spec = NULL,
                    gn_fit_spec = NULL,
                    Qn_fit_spec = NULL) {
  # check arguments and set up some objects for programmatic convenience
  call <- match.call(expand.dots = TRUE)
  estimator <- match.arg(estimator)
  fluctuation <- match.arg(fluctuation)
  eif_reg_type <- match.arg(eif_reg_type)

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

  # perform sub-setting of data and implement IPC weighting if required
  if (!all(C == 1) & !is.null(V)) {
    if (is.character(V)) {
      # combine censoring node information
      V_in <- data.table::as.data.table(mget(V))
      # NOTE: resolves downstream naming error
      V_names <- lapply(seq_along(V), function(j) {
        node <- mget(V[j], inherits = TRUE)[[1]]
        if (!is.null(dim(node))) {
          colnames(node)
        } else {
          V[j]
        }
      })
      colnames(V_in) <- do.call(c, V_names)
    } else {
      # assume V is a given matrix-type object with correctly set names
      V_in <- V
    }

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
      ipcw_estim <- ipcw_fit_spec
    } else {
      ipcw_estim <- do.call(est_ipcw, ipcw_estim_args)
    }

    # extract IPC weights for censoring case and normalize weights
    cens_weights <- ipcw_estim$ipc_weights

    # remove column corresponding to indicator for censoring
    data_internal <- data.table::data.table(W, A, C, Y)
    data_internal <- data_internal[C == 1, ] # NOTE: subset forces copying :(
    data_internal[, C := NULL]
  } else {
    # if no censoring, we can just use IPC weights that are identically 1
    V_in <- NULL
    cens_weights <- C
    ipcw_estim <- list(pi_mech = rep(1, length(C)), ipc_weights = C[C == 1])
    data_internal <- data.table::data.table(W, A, Y)
  }

  # initial estimate of the treatment mechanism (propensity score)
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
    if (g_fit_type == "hal") {
      # since fitting a GLM, can safely remove all args related to SL
      g_fit_args <- g_fit_args[!stringr::str_detect(names(g_fit_args), "sl")]

      # reshape args to a list suitable to be passed to do.call
      gn_estim_args <- unlist(
        list(gn_estim_in, std_args = list(g_fit_args)),
        recursive = FALSE
      )
    } else if (g_fit_type == "sl") {
      # if fitting SL, we can discard all the standard non-sl3 arguments
      g_fit_args <- g_fit_args[stringr::str_detect(names(g_fit_args), "sl")]

      # reshapes list of args to make passing to do.call possible
      gn_estim_args <- unlist(list(gn_estim_in, g_fit_args), recursive = FALSE)
    }

    # pass the relevant args for computing the propensity score
    gn_estim <- do.call(est_g, gn_estim_args)
  }

  # initial estimate of the outcome regression
  if (!is.null(Qn_fit_spec) & Q_fit_type == "fit_spec") {
    Qn_estim <- Qn_fit_spec
  } else {
    # generate and reshape args to pass to function for outcome regression
    Qn_estim_in <- list(
      Y = data_internal$Y,
      A = data_internal$A,
      W = data_internal[, W_names, with = FALSE],
      delta = delta,
      ipc_weights = cens_weights,
      fit_type = Q_fit_type
    )
    Qn_estim_args <- unlist(list(Qn_estim_in, Q_fit_args), recursive = FALSE)

    # invoke function to estimate outcome regression
    Qn_estim <- do.call(est_Q, Qn_estim_args)
  }

  # initial estimate of the auxiliary covariate
  Hn_estim <- est_Hn(gn = gn_estim)

  # compute targeted maximum likelihood estimator
  if (estimator == "tmle") {
    tmle_fit <- tmle_txshift(
      data_internal = data_internal,
      C = C,
      V = V_in,
      delta = delta,
      ipcw_estim = ipcw_estim,
      Qn_estim = Qn_estim,
      Hn_estim = Hn_estim,
      fluctuation = fluctuation,
      max_iter = max_iter,
      eif_reg_type = eif_reg_type,
      ipcw_fit_args = ipcw_fit_args,
      ipcw_efficiency = ipcw_efficiency
    )

    # return output object created by TML estimation routine
    tmle_fit$call <- call
    return(tmle_fit)

    # compute one-step (augmented inverse probability weighted) estimator
  } else if (estimator == "onestep") {
    onestep_fit <- onestep_txshift(
      data_internal = data_internal,
      C = C,
      V = V_in,
      delta = delta,
      ipcw_estim = ipcw_estim,
      Qn_estim = Qn_estim,
      Hn_estim = Hn_estim,
      eif_reg_type = eif_reg_type,
      ipcw_fit_args = ipcw_fit_args,
      ipcw_efficiency = ipcw_efficiency
    )

    # return output object created by AIPW estimation routine
    onestep_fit$call <- call
    return(onestep_fit)
  }
}
