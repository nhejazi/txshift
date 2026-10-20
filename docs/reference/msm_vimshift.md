# Working marginal structural model for causal effects of an intervention grid

Working marginal structural model for causal effects of an intervention
grid

## Usage

``` r
msm_vimshift(
  W,
  A,
  C_cens = rep(1, length(Y)),
  Y,
  C_samp = rep(1, length(Y)),
  V = NULL,
  delta_grid = seq(-0.5, 0.5, 0.5),
  msm_form = list(type = "linear", knot = NA),
  estimator = c("tmle", "onestep"),
  weighting = c("identity", "variance"),
  ci_level = 0.95,
  ci_type = c("marginal", "simultaneous"),
  ...
)
```

## Arguments

- W:

  A `matrix`, `data.frame`, or similar containing a set of baseline
  covariates.

- A:

  A `numeric` vector corresponding to a treatment variable. The
  parameter of interest is defined as a location shift of this quantity.

- C_cens:

  A `numeric` indicator for whether a given observation was subject to
  censoring by way of loss to follow-up. The default assumes no
  censoring due to loss to follow-up.

- Y:

  A `numeric` vector of the observed outcomes.

- C_samp:

  A `numeric` indicator for whether a given observation was subject to
  censoring by being omitted from the second-stage sample, used to
  compute an inverse probability of censoring weighted estimator in such
  cases. The default assumes no censoring due to two-phase sampling.

- V:

  The covariates that are used in determining the sampling procedure
  that gives rise to censoring. The default is `NULL` and corresponds to
  scenarios in which there is no censoring (in which case all values in
  the preceding argument `C` must be uniquely 1. To specify this, pass
  in a NAMED `list` identifying variables amongst W, A, Y that are
  thought to have played a role in defining the sampling/censoring
  mechanism (C).

- delta_grid:

  A `numeric` vector giving the individual values of the shift parameter
  used in computing each of the estimates.

- msm_form:

  A `list` specifying the type of working MSM to fit to summarize the
  counterfactual means. The `list` has two components: (1) `"type"`,
  which may be either "linear" or "piecewise", and (2) `"knot"`, which,
  if specified, must be a value in `delta_grid`. See examples for its
  use.

- estimator:

  The type of estimator to be fit, either `"tmle"` for targeted maximum
  likelihood estimation or `"onestep"` for a one-step augmented inverse
  probability weighted (AIPW) estimator.

- weighting:

  Whether to weight each parameter estimate by the inverse of its
  variance (in order to improve stability of the resultant MSM fit) or
  to simply weight all parameter estimates equally. The default is the
  option `"identity"`, weighting all estimates identically.

- ci_level:

  A `numeric` indicating the desired coverage level of the confidence
  interval to be computed.

- ci_type:

  Whether to construct a simultaneous confidence band covering all
  parameter estimates at once or marginal confidence intervals covering
  each parameter estimate separately. The default is to construct
  marginal confidence intervals for each parameter estimate rather than
  a simultaneous confidence band.

- ...:

  Additional arguments to be passed to
  [`txshift`](https://codex.nimahejazi.org/txshift/reference/txshift.md).

## Value

A `list` containing estimates of the individual counterfactual means
over a grid in the shift parameters (`delta_grid`), alongside the
estimate of a marginal structural model that summarizes a trend through
these counterfactual means.

## Details

Computes estimates of the counterfactual mean over a grid of shift
stochastic interventions and fits a working marginal structural model to
summarize the trend through the counterfactual means as a function of
the specified shift intervention. The working marginal structural model
may be linear in the shift parameter or piecewise linear with a single
knot point. Provides support for two weighting schemes, may be used with
either of the one-step or TML estimators, and also allows the
construction of marginal or simultaneous confidence intervals.

## Examples

``` r
if (require("sl3")) {
  n_obs <- 100
  W <- as.numeric(replicate(1, rbinom(n_obs, 1, 0.5)))
  A <- as.numeric(rnorm(n_obs, mean = 2 * W, sd = 1))
  Y <- rbinom(n_obs, 1, plogis(2 * A - W))
  msm <- msm_vimshift(
    W = W, A = A, Y = Y, estimator = "tmle",
    g_exp_fit_args = list(
      fit_type = "sl",
      sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
    ),
    Q_fit_args = list(
      fit_type = "glm",
      glm_formula = "Y ~ ."
    ),
    delta_grid = seq(-1, 1, 0.25)
  )

  # fit a linear spline with knot at 0
  n_obs <- 100
  W <- as.numeric(replicate(1, rbinom(n_obs, 1, 0.5)))
  A <- as.numeric(rnorm(n_obs, mean = 2 * W, sd = 1))
  Y <- rbinom(n_obs, 1, plogis(0.1 * A * (A >= 0) - 3 * A * (A < 0) - W))
  msm <- msm_vimshift(
    W = W, A = A, Y = Y, estimator = "tmle",
    g_exp_fit_args = list(
      fit_type = "sl",
      sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
    ),
    Q_fit_args = list(
      fit_type = "glm",
      glm_formula = "Y ~ ."
    ),
    delta_grid = seq(-1, 1, 0.25),
    msm_form = list(type = "piecewise", knot = 0)
  )
}
#> Loading required package: sl3
```
