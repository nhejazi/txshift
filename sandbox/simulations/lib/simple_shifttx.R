simple_shifttx_sim <- function(iter, samp, delta = 0.5,
                               n_w, w1_prob, a1_mean, a0_mean,
                               fit_type = "sl",
                               glm_form = "Y ~ .",
                               sl_lrnrs = c("mean", "glm_fast"),
                               sl_meta = "nnls") {

  if ((iter %% 100) == 0) {
    print(paste("Starting simulation", iter, "for sample size", samp))
  }

  ## baseline covariate -- simple, binary
  W <- as.numeric(replicate(n_w, rbinom(samp, 1, w1_prob)))

  ## set and organize treatment based on baseline W
  A1 <- rnorm(length(which(W == 1)), mean = a1_mean, sd = 1)
  A0 <- rnorm(length(which(W == 0)), mean = a0_mean, sd = 1)
  A <- rep(NA, samp)
  A[which(W == 0)] <- A0
  A[which(W == 1)] <- A1

  # create outcome
  Y <- A + W + rnorm(samp)

  # compute the truth from the simulation
  # \Psi(P_0) = E_0( \bar{Q}_0(A + \delta, W))
  #   = E_0(A + \delta + W)  -> definition of \bar{Q}_0
  #   = E_0(E_0(A \mid W) + \delta + W)  -> iterated expectation
  #   = \delta + P(W = 0) * (E_0(A \mid W = 0) + 0)
  #     + P(W = 1) * (E_0(A \mid W = 1) + 1)  -> defn of expectation over W
  #   = 0.5 + 0.5 * (0 + 0) + 0.5 * (2 + 1)
  #   = 2
  sim_truth = delta + (1 - w1_prob) * (a0_mean + 0) + w1_prob * (a1_mean + 1)

  # do a TMLE
  tmle_shift <- tmle_shifttx(W = W, A = A, Y = Y, delta = delta,
                             g_fit_args = list(nbins = 20,
                                               bin_method = "dhist",
                                               bin_estimator =
                                                 speedglmR6$new(),
                                               parfit = FALSE),
                             Q_fit_args = list(fit_method = fit_type,
                                               glm_formula = glm_form,
                                               sl_learners = sl_lrnrs,
                                               sl_metalearner = sl_meta),
                             fluc_method = "standard",
                             eif_tol = 1e-7
                            )
  ci_shift <- confint(tmle_shift)

  # summary statistics from TMLE procedure
  point_est <- tmle_shift$psi
  ci_in <- as.numeric(between(sim_truth, ci_shift[1], ci_shift[3]))
  sim_iter_out <- c(ci_shift[1], point_est, ci_shift[3], ci_in, sim_truth)
  names(sim_iter_out) <- c("lwr_ci", "est_psi", "upr_ci", "truth_in_ci",
                           "truth")
  sim_iter_out <- as.data.frame(t(as.matrix(sim_iter_out)))
  return(sim_iter_out)
}

