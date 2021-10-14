#' Efficient Estimate of Counterfactual Mean of Stochastic Shift Intervention
#'
#' @details Construct a one-step estimate or targeted minimum loss estimate of
#'  the counterfactual mean under a modified treatment policy, automatically
#'  making adjustments for two-phase sampling when a censoring indicator is
#'  included. Ensemble machine learning may be used to construct the initial
#'  estimates of nuisance functions using \pkg{sl3}.
#'
#' @param W A \code{matrix}, \code{data.frame}, or similar containing a set of
#'  baseline covariates.
#' @param A A \code{numeric} vector corresponding to a treatment variable. The
#'  parameter of interest is defined as a location shift of this quantity.
#' @param C_cens A \code{numeric} indicator for whether a given observation was
#'  subject to censoring by way of loss to follow-up. The default assumes no
#'  censoring due to loss to follow-up.
#' @param Y A \code{numeric} vector of the observed outcomes.
#' @param C_samp A \code{numeric} indicator for whether a given observation was
#'  subject to censoring by being omitted from the second-stage sample, used to
#'  compute an inverse probability of censoring weighted estimator in such
#'  cases. The default assumes no censoring due to two-phase sampling.
#' @param V The covariates that are used in determining the sampling procedure
#'  that gives rise to censoring. The default is \code{NULL} and corresponds to
#'  scenarios in which there is no censoring (in which case all values in the
#'  preceding argument \code{C_samp} must be uniquely 1). To specify this, pass
#'  in a \code{character} vector identifying variables amongst W, A, Y thought
#'  to have impacted the definition of the sampling mechanism (C_samp). This
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
#' @param samp_fit_args A \code{list} of arguments, all but one of which are
#'  passed to \code{\link{est_samp}}. For details, consult the documentation of
#'  \code{\link{est_samp}}. The first element (i.e., \code{fit_type}) is used
#'  to determine how this regression is fit: generalized linear model ("glm")
#'  or Super Learner ("sl"), and "external" a user-specified input of the form
#'  produced by \code{\link{est_samp}}.
#' @param g_exp_fit_args A \code{list} of arguments, all but one of which are
#'  passed to \code{\link{est_g_exp}}. For details, see the documentation of
#'  \code{\link{est_g_exp}}. The 1st element (i.e., \code{fit_type}) specifies
#'  how this regression is fit: \code{"hal"} to estimate conditional densities
#'  via the highly adaptive lasso (via \pkg{haldensify}), \code{"sl"} for
#'  \pkg{sl3} learners used to fit Super Learner ensembles to densities via
#'  \pkg{sl3}'s \code{Lrnr_haldensify} or similar, and \code{"external"} for
#'  user-specified input of the form produced by \code{\link{est_g_exp}}.
#' @param g_cens_fit_args A \code{list} of arguments, all but one of which are
#'  passed to \code{\link{est_g_cens}}. For details, see the documentation of
#'  \code{\link{est_g_cens}}. The 1st element (i.e., \code{fit_type}) specifies
#'  how this regression is fit: \code{"glm"} for a generalized linear model
#'  or \code{"sl"} for \pkg{sl3} learners used to fit a Super Learner ensemble
#'  for the censoring mechanism, and \code{"external"} for user-specified input
#'  of the form produced by \code{\link{est_g_cens}}.
#' @param Q_fit_args A \code{list} of arguments, all but one of which are
#'  passed to \code{\link{est_Q}}. For details, consult the documentation for
#'  \code{\link{est_Q}}. The first element (i.e., \code{fit_type}) is used to
#'  determine how this regression is fit: \code{"glm"} for a generalized linear
#'  model for the outcome mechanism, \code{"sl"} for \pkg{sl3} learners used
#'  to fit a Super Learner for the outcome mechanism, and \code{"external"}
#'  for user-specified input of the form produced by \code{\link{est_Q}}.
#' @param eif_reg_type Whether a flexible nonparametric function ought to be
#'  used in the dimension-reduced nuisance regression of the targeting step for
#'  the censored data case. By default, the method used is a nonparametric
#'  regression based on the Highly Adaptive Lasso (from \pkg{hal9001}). Set
#'  this to \code{"glm"} to instead use a simple linear regression model. In
#'  this step, the efficient influence function (EIF) is regressed against
#'  covariates contributing to the censoring mechanism (i.e., EIF ~ V | C = 1).
#' @param ipcw_efficiency Whether to use an augmented inverse probability of
#'  censoring weighted EIF estimating equation to ensure efficiency of the
#'  resultant estimate. The default is \code{TRUE}; the inefficient estimation
#'  procedure specified by \code{FALSE} is only supported for completeness.
#' @param samp_fit_ext The results of an external fitting procedure used to
#'  estimate the two-phase sampling mechanism, to be used in constructing the
#'  inverse probability of censoring weighted TML or one-step estimator. The
#'  input provided must match the output of \code{\link{est_samp}} exactly.
#' @param gn_exp_fit_ext The results of an external fitting procedure used to
#'  estimate the exposure mechanism (generalized propensity score), to be used
#'  in constructing the TML or one-step estimator. The input provided must
#'  match the output of \code{\link{est_g_exp}} exactly.
#' @param gn_cens_fit_ext The results of an external fitting procedure used to
#'  estimate the censoring mechanism (propensity score for missingness), to be
#'  used in constructing the TML or one-step estimator. The input provided must
#'  match the output of \code{\link{est_g_cens}} exactly.
#' @param Qn_fit_ext The results of an external fitting procedure used to
#'  estimate the outcome mechanism, to be used in constructing the TML or
#'  one-step estimator. The input provided must match the output of
#'  \code{\link{est_Q}} exactly; use of this argument is only recommended for
#'  power users.
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
#' @export
#'
#' @examples
#' set.seed(429153)
#' n_obs <- 100
#' W <- replicate(2, rbinom(n_obs, 1, 0.5))
#' A <- rnorm(n_obs, mean = 2 * W, sd = 1)
#' Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
#' C_samp <- rbinom(n_obs, 1, plogis(W + Y)) # two-phase sampling
#' C_cens <- rbinom(n_obs, 1, plogis(rowSums(W) + 0.5))
#'
#' # construct a TML estimate, ignoring censoring
#' tmle <- txshift(
#'   W = W, A = A, Y = Y, delta = 0.5,
#'   estimator = "onestep",
#'   g_exp_fit_args = list(
#'     fit_type = "hal",
#'     n_bins = 3,
#'     lambda_seq = exp(seq(-1, -10, length = 50))
#'   ),
#'   Q_fit_args = list(
#'     fit_type = "glm",
#'     glm_formula = "Y ~ ."
#'   )
#' )
#'
#' \dontrun{
#' # construct a TML estimate, accounting for censoring
#' tmle <- txshift(
#'   W = W, A = A, C_cens = C_cens, Y = Y, delta = 0.5,
#'   estimator = "onestep",
#'   g_exp_fit_args = list(
#'     fit_type = "hal",
#'     n_bins = 3,
#'     lambda_seq = exp(seq(-1, -10, length = 50))
#'   ),
#'   g_cens_fit_args = list(
#'     fit_type = "glm",
#'     glm_formula = "C_cens ~ ."
#'   ),
#'   Q_fit_args = list(
#'     fit_type = "glm",
#'     glm_formula = "Y ~ ."
#'   )
#' )
#'
#' # construct a TML estimate under two-phase sampling, ignoring censoring
#' ipcwtmle <- txshift(
#'   W = W, A = A, Y = Y, delta = 0.5,
#'   C_samp = C_samp, V = c("W", "Y"),
#'   estimator = "onestep", max_iter = 3,
#'   samp_fit_args = list(fit_type = "glm"),
#'   g_exp_fit_args = list(
#'     fit_type = "hal",
#'     n_bins = 3,
#'     lambda_seq = exp(seq(-1, -10, length = 50))
#'   ),
#'   Q_fit_args = list(
#'     fit_type = "glm",
#'     glm_formula = "Y ~ ."
#'   ),
#'   eif_reg_type = "glm"
#' )
#' }
#'
#' # construct a TML estimate acconting for two-phase sampling and censoring
#' ipcwtmle <- txshift(
#'   W = W, A = A, C_cens = C_cens, Y = Y, delta = 0.5,
#'   C_samp = C_samp, V = c("W", "Y"),
#'   estimator = "onestep", max_iter = 3,
#'   samp_fit_args = list(fit_type = "glm"),
#'   g_exp_fit_args = list(
#'     fit_type = "hal",
#'     n_bins = 3,
#'     lambda_seq = exp(seq(-1, -10, length = 50))
#'   ),
#'   g_cens_fit_args = list(
#'     fit_type = "glm",
#'     glm_formula = "C_cens ~ ."
#'   ),
#'   Q_fit_args = list(
#'     fit_type = "glm",
#'     glm_formula = "Y ~ ."
#'   ),
#'   eif_reg_type = "glm"
#' )
txshift <- function(W,
                    A,
                    C_cens = rep(1, length(A)),
                    Y,
                    C_samp = rep(1, length(Y)),
                    V = NULL,
                    delta = 0,
                    estimator = c("tmle", "onestep"),
                    fluctuation = c("standard", "weighted"),
                    max_iter = 10,
                    samp_fit_args = list(
                      fit_type = c("glm", "sl", "external"),
                      sl_learners = NULL
                    ),
                    g_exp_fit_args = list(
                      fit_type = c("hal", "sl", "external"),
                      lambda_seq = exp(seq(-1, -13, length = 300)),
                      sl_learners_density = NULL
                    ),
                    g_cens_fit_args = list(
                      fit_type = c("glm", "sl", "external"),
                      glm_formula = "C_cens ~ .",
                      sl_learners = NULL
                    ),
                    Q_fit_args = list(
                      fit_type = c("glm", "sl", "external"),
                      glm_formula = "Y ~ .",
                      sl_learners = NULL
                    ),
                    eif_reg_type = c("hal", "glm"),
                    ipcw_efficiency = TRUE,
                    samp_fit_ext = NULL,
                    gn_exp_fit_ext = NULL,
                    gn_cens_fit_ext = NULL,
                    Qn_fit_ext = NULL) {
  # check arguments and set up some objects for programmatic convenience
  call <- match.call(expand.dots = TRUE)
  estimator <- match.arg(estimator)
  fluctuation <- match.arg(fluctuation)
  eif_reg_type <- match.arg(eif_reg_type)

  # dissociate fit type from other arguments to simplify passing to do.call
  samp_fit_type <- unlist(samp_fit_args[names(samp_fit_args) == "fit_type"],
    use.names = FALSE
  )
  g_exp_fit_type <- unlist(g_exp_fit_args[names(g_exp_fit_args) == "fit_type"],
    use.names = FALSE
  )
  g_cens_fit_type <-
    unlist(g_cens_fit_args[names(g_cens_fit_args) == "fit_type"],
      use.names = FALSE
    )
  Q_fit_type <- unlist(Q_fit_args[names(Q_fit_args) == "fit_type"],
    use.names = FALSE
  )
  samp_fit_args <- samp_fit_args[names(samp_fit_args) != "fit_type"]
  g_exp_fit_args <- g_exp_fit_args[names(g_exp_fit_args) != "fit_type"]
  g_cens_fit_args <- g_cens_fit_args[names(g_cens_fit_args) != "fit_type"]
  Q_fit_args <- Q_fit_args[names(Q_fit_args) != "fit_type"]

  # coerce W to matrix and, if no names in W, assign them generically
  if (!is.matrix(W)) W <- as.matrix(W)
  W_names <- colnames(W)
  if (is.null(W_names)) {
    W_names <- paste0("W", seq_len(ncol(W)))
    colnames(W) <- W_names
  }

  # subset data and implement IPC weighting for two-phase sampling corrections
  if (!all(C_samp == 1) && !is.null(V)) {
    if (is.character(V)) {
      # combine censoring node information
      V_in <- data.table::as.data.table(do.call(cbind, mget(V)))
    } else {
      # assume V is a given matrix-type object with correctly set names
      V_in <- V
    }

    # create arguments for sampling mechanism estimation
    samp_estim_in <- list(
      V = V_in, C_samp = C_samp,
      fit_type = samp_fit_type
    )

    # reshapes the list of args so that it can be passed to do.call
    samp_estim_args <- unlist(
      list(samp_estim_in, samp_fit_args),
      recursive = FALSE
    )

    # compute the IPC weights by passing all args to the relevant function
    if (!is.null(samp_fit_ext) && samp_fit_type == "external") {
      samp_estim <- samp_fit_ext
    } else {
      samp_estim <- do.call(est_samp, samp_estim_args)
    }

    # extract IPC weights for two-phase sampling
    samp_weights <- (C_samp / samp_estim)[C_samp == 1]

    # remove column corresponding to indicator for censoring
    data_internal <- data.table::data.table(W, A, C_cens, Y, C_samp)
    data_internal <- data_internal[C_samp == 1, ] # subset forced copy :(
    data_internal[, C_samp := NULL]
  } else {
    # if no two-phase sampling, we can use IPC weights that are identically 1
    V_in <- NULL
    samp_weights <- samp_estim <- C_samp
    data_internal <- data.table::data.table(W, A, C_cens, Y)
  }

  # initial estimate of the treatment mechanism (generalized propensity score)
  if (!is.null(gn_exp_fit_ext) && g_exp_fit_type == "external") {
    gn_exp_estim <- gn_exp_fit_ext
  } else {
    gn_exp_estim_in <- list(
      A = data_internal$A,
      W = data_internal[, W_names, with = FALSE],
      delta = delta,
      samp_weights = samp_weights,
      fit_type = g_exp_fit_type
    )
    if (g_exp_fit_type == "hal") {
      # reshape args to a list suitable to be passed to do.call
      gn_exp_estim_args <- unlist(
        list(gn_exp_estim_in, haldensify_args = list(g_exp_fit_args)),
        recursive = FALSE
      )
    } else if (g_exp_fit_type == "sl") {
      # reshapes list of args to make passing to do.call possible
      gn_exp_estim_args <- unlist(list(gn_exp_estim_in, g_exp_fit_args),
        recursive = FALSE
      )
    }

    # pass the relevant args for computing the propensity score
    gn_exp_estim <- do.call(est_g_exp, gn_exp_estim_args)
  }

  # estimate the natural censoring mechanism for joint intervention
  if (any(C_cens != 1)) {
    if (!is.null(gn_cens_fit_ext) && g_cens_fit_type == "external") {
      gn_cens_weights <- C_cens[C_samp == 1] / gn_cens_fit_ext
    } else {
      gn_cens_estim_in <- list(
        C = data_internal$C_cens,
        A = data_internal$A,
        W = data_internal[, W_names, with = FALSE],
        samp_weights = samp_weights,
        fit_type = g_cens_fit_type
      )
      gn_cens_estim_args <- unlist(list(gn_cens_estim_in, g_cens_fit_args),
        recursive = FALSE
      )

      # invoke function to estimate the natural censoring mechanism
      gn_cens_estim <- do.call(est_g_cens, gn_cens_estim_args)
      gn_cens_weights <- C_cens[C_samp == 1] / gn_cens_estim
    }
  } else {
    gn_cens_weights <- rep(1, nrow(data_internal))
  }

  # initial estimate of the outcome mechanism
  if (!is.null(Qn_fit_ext) && Q_fit_type == "external") {
    Qn_estim <- Qn_fit_ext
  } else {
    # generate and reshape args to pass to function for outcome regression
    Qn_estim_in <- list(
      Y = data_internal$Y,
      C_cens = data_internal$C_cens,
      A = data_internal$A,
      W = data_internal[, W_names, with = FALSE],
      delta = delta,
      samp_weights = samp_weights,
      fit_type = Q_fit_type
    )
    Qn_estim_args <- unlist(list(Qn_estim_in, Q_fit_args), recursive = FALSE)

    # invoke function to estimate outcome regression
    Qn_estim <- do.call(est_Q, Qn_estim_args)
  }

  # initial estimate of the auxiliary covariate
  Hn_estim <- est_Hn(gn_exp = gn_exp_estim)

  # compute whichever efficient estimator was asked for
  if (estimator == "tmle") {
    # compute targeted maximum likelihood estimator
    tmle_fit <- tmle_txshift(
      data_internal = data_internal,
      C_samp = C_samp,
      V = V_in,
      delta = delta,
      samp_estim = samp_estim,
      gn_cens_weights = gn_cens_weights,
      Qn_estim = Qn_estim,
      Hn_estim = Hn_estim,
      fluctuation = fluctuation,
      max_iter = max_iter,
      eif_reg_type = eif_reg_type,
      samp_fit_args = samp_fit_args,
      ipcw_efficiency = ipcw_efficiency
    )

    # return output object created by the TML estimation routine
    tmle_fit$call <- call
    return(tmle_fit)
  } else if (estimator == "onestep") {
    # compute the efficient one-step estimator
    onestep_fit <- onestep_txshift(
      data_internal = data_internal,
      C_samp = C_samp,
      V = V_in,
      delta = delta,
      samp_estim = samp_estim,
      gn_cens_weights = gn_cens_weights,
      Qn_estim = Qn_estim,
      Hn_estim = Hn_estim,
      eif_reg_type = eif_reg_type,
      samp_fit_args = samp_fit_args,
      ipcw_efficiency = ipcw_efficiency
    )

    # return output object created by the one-step estimation routine
    onestep_fit$call <- call
    return(onestep_fit)
  }
}
