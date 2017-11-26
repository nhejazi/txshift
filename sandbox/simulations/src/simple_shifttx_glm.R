# infrastructure packages and setup
library(ProjectTemplate)
load.project()
library(future)
library(doFuture)
registerDoFuture()
library(data.table)
library(here)
source(here("lib", "simple_shifttx.R"))

# estimation packages
library(sl3)
library(condensier)
library(shifttx)

# meta-parameters for controlling simulation type
seed_int <- 9461287
n_obs <- c(50, 250, 1000)
n_sim <- 1e3


# parameters for simulating simple data for tmle-shift sketch
n_w <- 1
w_win_prob <- 0.5
a1_mean <- 2
a0_mean <- 0
delta_shift <- 0.5


### Fitting outcome model with GLMs

# simulation
simple_sim_glm <- foreach(samp_it = seq_along(n_obs), .combine = rbind) %do% {
  # set sample size from foreach loop
  samp <- n_obs[samp_it]

  # run simulation for given sample size
  sim_results <- foreach(i = seq_len(n_sim),
                         .combine = rbind) %dopar% {
    # set seed in foreach
    set.seed(seed_int + i)

    # run the TMLE procedure
    out <- simple_shifttx_sim(iter = i, samp = samp, delta = delta_shift,
                              n_w = n_w, w1_prob = w_win_prob,
                              a1_mean = a1_mean, a0_mean = a0_mean,
                              fit_type = "glm",
                              glm_form = "Y ~ .")
    out
  }
  # compute average estimate, bias, and sd across simulations
  sim_point_est <- mean(sim_results$est_psi)
  sim_est_var <- var(sim_results$est_psi)
  sim_est_bias <- sim_point_est - mean(sim_results$truth)
  sim_ci_cover <- sum(sim_results$truth_in_ci) / nrow(sim_results)
  sim_out <- c(samp, sim_point_est, sim_est_var, sim_est_bias, sim_ci_cover)
  names(sim_out) <- c("n_size", "est", "var", "bias", "coverage")
  sim_samp_out <- as.data.frame(t(as.matrix(sim_out)))
  sim_samp_out
}

simple_sim_glm

