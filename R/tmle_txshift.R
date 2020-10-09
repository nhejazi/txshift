#' Targeted Minimum Loss Estimate of Counterfactual Mean of Stochastic Shift
#' Intervention
#'
#' @details Invokes the procedure to construct a targeted minimum loss estimate
#'  (TMLE) of the counterfactual mean under a modified treatment policy.
#'
#' @param data_internal A \code{data.table} constructed internally by a call to
#'  \code{\link{txshift}}. This contains most of the data for computing the
#'  targeted minimum loss (TML) estimator.
#' @param C_samp A \code{numeric} indicator for whether a given observation was
#'  included in the second-stage sample, used to compute an IPC-weighted
#'  one-step estimator in cases where two-stage sampling is performed. Default
#'  assumes no censoring due to sampling.
#' @param V The covariates that are used in determining the sampling procedure
#'  that gives rise to censoring. The default is \code{NULL} and corresponds to
#'  scenarios in which there is no censoring (in which case all values in the
#'  preceding argument \code{C_samp} must be 1. To specify this, pass in a
#'  NAMED \code{list} identifying variables amongst W, A, Y that are thought to
#'  have played a role in defining the sampling mechanism.
#' @param delta A \code{numeric} value indicating the shift in the treatment to
#'  be used in defining the target parameter. This is defined with respect to
#'  the scale of the treatment (A).
#' @param samp_estim An object providing the value of the sampling mechanism
#'  evaluated across the full data. This object is passed in after being
#'  constructed by a call to the internal function \code{\link{est_samp}}.
#' @param gn_cens_weights TODO: document
#' @param Qn_estim An object providing the value of the outcome evaluated after
#'  imposing a shift in the treatment. This object is passed in after being
#'  constructed by a call to the internal function \code{\link{est_Q}}.
#' @param Hn_estim An object providing values of the auxiliary ("clever")
#'  covariate, constructed from the treatment mechanism and required for
#'  targeted minimum loss-based estimation. This object object should be passed
#'  in after being constructed by a call to \code{\link{est_Hn}}.
#' @param fluctuation The method to be used in the submodel fluctuation step
#'  (targeting step) to compute the TML estimator. The choices are "standard"
#'  and "weighted" for where to place the auxiliary covariate in the logistic
#'  tilting regression.
#' @param max_iter A \code{numeric} integer giving the maximum number of steps
#'  to be taken in iterating to a solution of the efficient influence function.
#' @param eif_reg_type Whether a flexible nonparametric function ought to be
#'  used in the dimension-reduced nuisance regression of the targeting step for
#'  the censored data case. By default, the method used is a nonparametric
#'  regression based on the Highly Adaptive Lasso (from \pkg{hal9001}).
#'  Set this to \code{"glm"} to instead use a simple linear regression model.
#'  In this step, the efficient influence function (EIF) is regressed against
#'  covariates contributing to the censoring mechanism (i.e., EIF ~ V | C = 1).
#' @param samp_fit_args A \code{list} of arguments, all but one of which are
#'  passed to \code{\link{est_samp}}. For details, consult the documentation
#'  for \code{\link{est_samp}}. The first element (i.e., \code{fit_type}) is
#'  used to determine how this regression is fit: "glm" for generalized linear
#'  model, "sl" for a Super Learner, and "external" for a user-specified input
#'  of the form produced by \code{\link{est_samp}}.
#' @param ipcw_efficiency Whether to invoke an augmentation of the IPCW-TMLE
#'  procedure that performs an iterative process to ensure efficiency of the
#'  resulting estimate. The default is \code{TRUE}; set to \code{FALSE} to use
#'  an IPC-weighted loss rather than the IPC-augmented influence function.
#'
#' @importFrom assertthat assert_that
#' @importFrom data.table as.data.table setnames
#' @importFrom stringr str_detect
#' @importFrom Rdpack reprompt
#'
#' @return S3 object of class \code{txshift} containing the results of the
#'  procedure to compute a TML estimate of the treatment shift parameter.
tmle_txshift <- function(data_internal,
                         C_samp = rep(1, nrow(data_internal)),
                         V = NULL,
                         delta = 0,
                         samp_estim,
                         gn_cens_weights,
                         Qn_estim,
                         Hn_estim,
                         fluctuation = c("standard", "weighted"),
                         max_iter = 10,
                         eif_reg_type = c("hal", "glm"),
                         samp_fit_args,
                         ipcw_efficiency = TRUE) {
  # initialize counter
  n_steps <- 0

  # normalize (two-phase) sampling mechanism weights
  samp_weights <- C_samp / samp_estim
  samp_weights_norm <- samp_weights / sum(samp_weights)

  # invoke efficient IPCW-TMLE if satisfied; otherwise ineffecient variant
  if (!all(C_samp == 1) && ipcw_efficiency) {
    browser()
    # checks for necessary components for augmented EIF procedure
    assertthat::assert_that(!is.null(samp_estim))
    assertthat::assert_that(!is.null(V))

    # programmatic bookkeeping
    eif_mean <- Inf
    eif_tol <- 1 / length(samp_weights)
    conv_res <- matrix(replicate(3, rep(NA_real_, max_iter)), nrow = max_iter)

    # quantities to be updated across targeting iterations
    pi_mech_star <- samp_estim
    Qn_estim_updated <- Qn_estim

    # iterate procedure until convergence conditions are satisfied
    while (abs(eif_mean) > eif_tol && n_steps <= max_iter) {
      # iterate counter
      n_steps <- n_steps + 1

      # update fluctuation model, re-compute EIF, overwrite EIF
      tmle_ipcw_eif <- ipcw_eif_update(
        data_in = data_internal,
        C = C,
        V = V,
        ipc_mech = pi_mech_star,
        ipc_weights = samp_weights,
        ipc_weights_norm = samp_weights_norm,
        Qn_estim = Qn_estim_updated,
        Hn_estim = Hn_estim, # N.B., g_n never gets updated in this procedure
        estimator = "tmle",
        fluctuation = fluctuation,
        eif_reg_type = eif_reg_type
      )

      # overwrite and update quantities to be used in the next iteration
      Qn_estim_updated <- data.table::as.data.table(
        list(
          # NOTE: need to re-scale estimated Q_n within bounds of Y
          scale_to_unit(
            vals = tmle_ipcw_eif$fluc_mod_out$Qn_noshift_star
          ),
          scale_to_unit(
            vals = tmle_ipcw_eif$fluc_mod_out$Qn_shift_star
          )
        )
      )
      data.table::setnames(Qn_estim_updated, names(Qn_estim))

      # updated sampling/censoring weights and stabilize
      samp_weights <- tmle_ipcw_eif$ipc_weights
      samp_weights <- samp_weights / mean(samp_weights)
      samp_weights_norm <- tmle_ipcw_eif$ipc_weights_norm
      pi_mech_star <- tmle_ipcw_eif$pi_mech_star

      # compute updated mean of efficient influence function and save
      eif_ipcw <- tmle_ipcw_eif$eif_eval$eif - tmle_ipcw_eif$ipcw_eif_component
      eif_mean <- mean(eif_ipcw)
      eif_var <- var(eif_ipcw) / length(samp_weights)
      conv_res[n_steps, ] <- c(tmle_ipcw_eif$eif_eval$psi, eif_var, eif_mean)

      # TMLE convergence criterion based on re-scaled standard error
      tol_scaling <- 1 / (max(10, log(length(samp_weights))))
      eif_tol <- sqrt(eif_var) * tol_scaling
    }
    conv_res <- data.table::as.data.table(conv_res)
    data.table::setnames(conv_res, c("psi", "var", "eif_mean"))

    # return only the useful convergence results
    conv_res <- conv_res[!is.na(rowSums(conv_res)), ]

    # create output object
    txshift_out <- unlist(
      list(
        psi = as.numeric(unlist(conv_res[n_steps, "psi"])),
        var = as.numeric(unlist(conv_res[n_steps, "var"])),
        eif = list(eif_ipcw),
        iter_res = list(conv_res),
        n_iter = n_steps,
        estimator = "tmle",
        outcome = list(data_internal$Y)
      ),
      recursive = FALSE
    )
  } else {
    # fit logistic regression to fluctuate along the sub-model
    fitted_fluc_mod <- fit_fluctuation(
      Y = data_internal$Y,
      Qn_scaled = Qn_estim,
      Hn = Hn_estim,
      ipc_weights = (gn_cens_weights * samp_weights),
      method = fluctuation
    )

    # compute TML estimate and EIF for the treatment shift parameter
    tmle_eif_out <- eif(
      Y = data_internal$Y,
      Qn = Qn_estim,
      Hn = Hn_estim,
      estimator = "tmle",
      fluc_mod_out = fitted_fluc_mod,
      C_samp = C_samp,
      ipc_weights = (gn_cens_weights * samp_weights)
    )

    # create output object
    txshift_out <- unlist(
      list(
        psi = tmle_eif_out[["psi"]],
        var = tmle_eif_out[["var"]],
        eif = list(tmle_eif_out[["eif"]]),
        iter_res = NULL,
        n_iter = n_steps,
        estimator = "tmle",
        outcome = list(data_internal$Y)
      ),
      recursive = FALSE
    )
  }

  # S3-ify and return output object
  class(txshift_out) <- "txshift"
  return(txshift_out)
}

###############################################################################

#' Fit One-Dimensional Fluctuation Model for Updating Initial Estimates
#'
#' @details Procedure for fitting a one-dimensional fluctuation model to update
#'  the initial estimates of the outcome regression based on the auxiliary
#'  covariate. These updated estimates are subsequently used to construct the
#'  TML estimator of the counterfactual mean under a modified treatment policy.
#'
#' @param Y A \code{numeric} vector corresponding to an outcome variable.
#' @param Qn_scaled An object providing the value of the outcome evaluate
#'  after inducing a shift in the exposure. This object should be passed in
#'  after being constructed by a call to \code{\link{est_Q}}.
#' @param Hn An object providing values of the auxiliary ("clever") covariate,
#'  constructed from the treatment mechanism and required for targeted minimum
#'  loss estimation. This object object should be passed in after being
#'  constructed by a call to \code{\link{est_Hn}}.
#' @param ipc_weights A \code{numeric} vector that gives inverse probability of
#'  censoring weights for each observation. These are generated by invoking the
#'  routines for estimating the censoring mechanism.
#' @param method A \code{character} giving the type of regression to be used in
#'  traversing the fluctuation sub-model. The available choices are "weighted"
#'  and "standard". Consult the literature for details on the differences.
#' @param flucmod_tol A \code{numeric} indicating the largest value to be
#'  tolerated in the fluctuation model for the targeted minimum loss estimator.
#'
#' @importFrom stats qlogis glm fitted predict as.formula coef
#' @importFrom data.table as.data.table setnames
#'
#' @return A \code{list} containing the fluctuation model (a \code{glm} object)
#'  produced by logistic regression, a \code{character} vector indicating the
#'  type of fluctuation (whether the auxiliary covariates was used as a weight
#'  or included directly in the model formula), the updated estimates of the
#'  outcome regression under the shifted value of the exposure, and the updated
#'  estimates of the outcome regression under the natural value of exposure.
fit_fluctuation <- function(Y,
                            Qn_scaled,
                            Hn,
                            ipc_weights = rep(1, length(Y)),
                            method = c("standard", "weighted"),
                            flucmod_tol = 50) {

  # scale the outcome for the logit transform
  y_star <- scale_to_unit(
    vals = Y
  )

  # bound precision for use of logit transform
  Qn_scaled_bounded <- data.table::as.data.table(apply(
    Qn_scaled, 2,
    bound_precision
  ))

  # extract Q and obtain logit transform
  Qn_noshift_logit <- stats::qlogis(Qn_scaled_bounded$noshift)

  # fit the fluctuation regression in one of two ways
  if (method == "standard") {
    # note that \epsilon_n will be the coefficient of the covariate Hn
    suppressWarnings(
      fluctuation_model <- stats::glm(
        formula = stats::as.formula(
          "y_star ~ -1 + offset(logit_Qn) + Hn"
        ),
        data = data.table::as.data.table(list(
          y_star = y_star,
          logit_Qn = Qn_noshift_logit,
          Hn = Hn$noshift
        )),
        weights = ipc_weights,
        family = "binomial",
        start = 0
      )
    )
    coefs_fluc <- stats::coef(fluctuation_model)

    # check convergence of fluctuation model and sanity of estimates
    if (!fluctuation_model$converged || abs(max(coefs_fluc)) > flucmod_tol) {
      suppressWarnings(
        fluctuation_model <- stats::glm(
          formula = stats::as.formula("y_star ~ -1 + offset(logit_Qn) + Hn"),
          data = data.table::as.data.table(list(
            y_star = y_star,
            logit_Qn = Qn_noshift_logit,
            Hn = Hn$noshift
          )),
          weights = ipc_weights,
          family = "binomial"
        )
      )
      # if the fluctuation model hasn't converged or is unstable, simply set
      # the coefficients to disable updating, i.e., coef(Hn) := 0
      if (!fluctuation_model$converged | abs(max(coefs_fluc)) > flucmod_tol) {
        fluctuation_model$coefficients <- 0
      }
    }
  } else if (method == "weighted") {
    # note that epsilon_n will be the intercept term here
    suppressWarnings(
      fluctuation_model <- stats::glm(
        formula = stats::as.formula("y_star ~ offset(logit_Qn)"),
        data = data.table::as.data.table(list(
          y_star = y_star,
          logit_Qn = Qn_noshift_logit
        )),
        weights = as.numeric(Hn$noshift * ipc_weights),
        family = "binomial",
        start = 0
      )
    )
    coefs_fluc <- stats::coef(fluctuation_model)

    # check covergence of fluctuation model and sanity of estimates
    if (!fluctuation_model$converged || abs(max(coefs_fluc)) > flucmod_tol) {
      suppressWarnings(
        fluctuation_model <- stats::glm(
          formula = stats::as.formula("y_star ~ offset(logit_Qn)"),
          data = data.table::as.data.table(list(
            y_star = y_star,
            logit_Qn = Qn_noshift_logit
          )),
          weights = as.numeric(Hn$noshift * ipc_weights),
          family = "binomial"
        )
      )
      # if the updated fluctuation model hasn't converged or is unstable,
      # simply set the coefficient to zero to disable updating
      if (!fluctuation_model$converged | abs(max(coefs_fluc)) > flucmod_tol) {
        fluctuation_model$coefficients <- 0
      }
    }
  }

  # get fitted values from fluctuation model
  Qn_noshift_star_unit <- unname(stats::fitted(fluctuation_model))
  Qn_noshift_star <- scale_to_original(
    scaled_vals = Qn_noshift_star_unit,
    max_orig = max(Y),
    min_orig = min(Y)
  )

  # need to logit transform Qn(d(A,W),W)
  Qn_shift_logit <- stats::qlogis(Qn_scaled_bounded$upshift)

  # get Qn_star for the SHIFTED data
  if (method == "standard") {
    Qn_shift_star_data <- data.table::as.data.table(list(
      logit_Qn = Qn_shift_logit,
      Hn = Hn$shift
    ))

    # predict from fluctuation model on Q(d(A,W),W) and Hn(d(A,W))
    Qn_shift_star_unit <- unname(stats::predict(
      object = fluctuation_model,
      newdata = Qn_shift_star_data,
      type = "response"
    ))
  } else if (method == "weighted") {
    Qn_shift_star_data <- data.table::as.data.table(Qn_shift_logit)
    data.table::setnames(Qn_shift_star_data, "logit_Qn")

    # predict from fluctuation model on Q(d(A,W),W)
    Qn_shift_star_unit <- unname(stats::predict(
      object = fluctuation_model,
      newdata = Qn_shift_star_data,
      type = "response"
    ))
  }

  Qn_shift_star <- scale_to_original(
    scaled_vals = Qn_shift_star_unit,
    max_orig = max(Y),
    min_orig = min(Y)
  )

  # return the fit model object
  out <- list(
    fluc_fit = fluctuation_model,
    covar_method = method,
    Qn_shift_star = as.numeric(Qn_shift_star),
    Qn_noshift_star = as.numeric(Qn_noshift_star)
  )
  return(out)
}
