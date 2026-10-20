# Efficient Estimate of Counterfactual Mean of Stochastic Shift Intervention

Efficient Estimate of Counterfactual Mean of Stochastic Shift
Intervention

## Usage

``` r
txshift(
  W,
  A,
  C_cens = rep(1, length(A)),
  Y,
  C_samp = rep(1, length(Y)),
  V = NULL,
  delta = 0,
  estimator = c("tmle", "onestep"),
  fluctuation = c("standard", "weighted"),
  max_iter = 10,
  gps_bound = 0.01,
  samp_fit_args = list(fit_type = c("glm", "sl", "external"), sl_learners = NULL),
  g_exp_fit_args = list(fit_type = c("hal", "sl", "external"), lambda_seq = exp(seq(-1,
    -13, length = 300)), sl_learners_density = NULL),
  g_cens_fit_args = list(fit_type = c("glm", "sl", "external"), glm_formula =
    "C_cens ~ .^2", sl_learners = NULL),
  Q_fit_args = list(fit_type = c("glm", "sl", "external"), glm_formula = "Y ~ .^2",
    sl_learners = NULL),
  eif_reg_type = c("hal", "glm"),
  ipcw_efficiency = TRUE,
  samp_fit_ext = NULL,
  gn_exp_fit_ext = NULL,
  gn_cens_fit_ext = NULL,
  Qn_fit_ext = NULL
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
  the preceding argument `C_samp` must be uniquely 1). To specify this,
  pass in a `character` vector identifying variables amongst W, A, Y
  thought to have impacted the definition of the sampling mechanism
  (C_samp). This argument also accepts a `data.table` (or similar)
  object composed of combinations of variables W, A, Y; use of this
  option is NOT recommended.

- delta:

  A `numeric` value indicating the shift in the treatment to be used in
  defining the target parameter. This is defined with respect to the
  scale of the treatment (A).

- estimator:

  The type of estimator to be fit, either `"tmle"` for targeted maximum
  likelihood or `"onestep"` for a one-step estimator.

- fluctuation:

  The method to be used in the submodel fluctuation step (targeting
  step) to compute the TML estimator. The choices are "standard" and
  "weighted" for where to place the auxiliary covariate in the logistic
  tilting regression.

- max_iter:

  A `numeric` integer giving the maximum number of steps to be taken in
  iterating to a solution of the efficient influence function.

- gps_bound:

  `numeric` giving the lower limit of the generalized propensity score
  estimates to be tolerated (default = 0.01). Estimates falling below
  this value are truncated to this or 1/n. For details, see
  [`bound_propensity`](https://codex.nimahejazi.org/txshift/reference/bound_propensity.md).

- samp_fit_args:

  A `list` of arguments, all but one of which are passed to
  [`est_samp`](https://codex.nimahejazi.org/txshift/reference/est_samp.md).
  For details, consult the documentation of
  [`est_samp`](https://codex.nimahejazi.org/txshift/reference/est_samp.md).
  The first element (i.e., `fit_type`) is used to determine how this
  regression is fit: generalized linear model ("glm") or Super Learner
  ("sl"), and "external" a user-specified input of the form produced by
  [`est_samp`](https://codex.nimahejazi.org/txshift/reference/est_samp.md).

- g_exp_fit_args:

  A `list` of arguments, all but one of which are passed to
  [`est_g_exp`](https://codex.nimahejazi.org/txshift/reference/est_g_exp.md).
  For details, see the documentation of
  [`est_g_exp`](https://codex.nimahejazi.org/txshift/reference/est_g_exp.md).
  The 1st element (i.e., `fit_type`) specifies how this regression is
  fit: `"hal"` to estimate conditional densities via the highly adaptive
  lasso (via haldensify), `"sl"` for sl3 learners used to fit Super
  Learner ensembles to densities via sl3's `Lrnr_haldensify` or similar,
  and `"external"` for user-specified input of the form produced by
  [`est_g_exp`](https://codex.nimahejazi.org/txshift/reference/est_g_exp.md).

- g_cens_fit_args:

  A `list` of arguments, all but one of which are passed to
  [`est_g_cens`](https://codex.nimahejazi.org/txshift/reference/est_g_cens.md).
  For details, see the documentation of
  [`est_g_cens`](https://codex.nimahejazi.org/txshift/reference/est_g_cens.md).
  The 1st element (i.e., `fit_type`) specifies how this regression is
  fit: `"glm"` for a generalized linear model or `"sl"` for sl3 learners
  used to fit a Super Learner ensemble for the censoring mechanism, and
  `"external"` for user-specified input of the form produced by
  [`est_g_cens`](https://codex.nimahejazi.org/txshift/reference/est_g_cens.md).

- Q_fit_args:

  A `list` of arguments, all but one of which are passed to
  [`est_Q`](https://codex.nimahejazi.org/txshift/reference/est_Q.md).
  For details, consult the documentation for
  [`est_Q`](https://codex.nimahejazi.org/txshift/reference/est_Q.md).
  The first element (i.e., `fit_type`) is used to determine how this
  regression is fit: `"glm"` for a generalized linear model for the
  outcome mechanism, `"sl"` for sl3 learners used to fit a Super Learner
  for the outcome mechanism, and `"external"` for user-specified input
  of the form produced by
  [`est_Q`](https://codex.nimahejazi.org/txshift/reference/est_Q.md).

- eif_reg_type:

  Whether a flexible nonparametric function ought to be used in the
  dimension-reduced nuisance regression of the targeting step for the
  censored data case. By default, the method used is a nonparametric
  regression based on the Highly Adaptive Lasso (from hal9001). Set this
  to `"glm"` to instead use a simple linear regression model. In this
  step, the efficient influence function (EIF) is regressed against
  covariates contributing to the censoring mechanism (i.e., EIF ~ V \| C
  = 1).

- ipcw_efficiency:

  Whether to use an augmented inverse probability of censoring weighted
  EIF estimating equation to ensure efficiency of the resultant
  estimate. The default is `TRUE`; the inefficient estimation procedure
  specified by `FALSE` is only supported for completeness.

- samp_fit_ext:

  The results of an external fitting procedure used to estimate the
  two-phase sampling mechanism, to be used in constructing the inverse
  probability of censoring weighted TML or one-step estimator. The input
  provided must match the output of
  [`est_samp`](https://codex.nimahejazi.org/txshift/reference/est_samp.md)
  exactly.

- gn_exp_fit_ext:

  The results of an external fitting procedure used to estimate the
  exposure mechanism (generalized propensity score), to be used in
  constructing the TML or one-step estimator. The input provided must
  match the output of
  [`est_g_exp`](https://codex.nimahejazi.org/txshift/reference/est_g_exp.md)
  exactly.

- gn_cens_fit_ext:

  The results of an external fitting procedure used to estimate the
  censoring mechanism (propensity score for missingness), to be used in
  constructing the TML or one-step estimator. The input provided must
  match the output of
  [`est_g_cens`](https://codex.nimahejazi.org/txshift/reference/est_g_cens.md)
  exactly.

- Qn_fit_ext:

  The results of an external fitting procedure used to estimate the
  outcome mechanism, to be used in constructing the TML or one-step
  estimator. The input provided must match the output of
  [`est_Q`](https://codex.nimahejazi.org/txshift/reference/est_Q.md)
  exactly; use of this argument is only recommended for power users.

## Value

S3 object of class `txshift` containing the results of the procedure to
compute a TML or one-step estimate of the counterfactual mean under a
modified treatment policy that shifts a continuous-valued exposure by a
scalar amount `delta`. These estimates can be augmented to be consistent
and efficient when two-phase sampling is performed.

## Details

Construct a one-step estimate or targeted minimum loss estimate of the
counterfactual mean under a modified treatment policy, automatically
making adjustments for two-phase sampling when a censoring indicator is
included. Ensemble machine learning may be used to construct the initial
estimates of nuisance functions using sl3.

## Examples

``` r
set.seed(429153)
n_obs <- 100
W <- replicate(2, rbinom(n_obs, 1, 0.5))
A <- rnorm(n_obs, mean = 2 * W, sd = 1)
Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
C_samp <- rbinom(n_obs, 1, plogis(W + Y)) # two-phase sampling
C_cens <- rbinom(n_obs, 1, plogis(rowSums(W) + 0.5))

# construct a TML estimate, ignoring censoring
tmle <- txshift(
  W = W, A = A, Y = Y, delta = 0.5,
  estimator = "onestep",
  g_exp_fit_args = list(
    fit_type = "hal",
    n_bins = 3,
    lambda_seq = exp(seq(-1, -10, length = 50))
  ),
  Q_fit_args = list(
    fit_type = "glm",
    glm_formula = "Y ~ ."
  )
)
if (FALSE) { # \dontrun{
# construct a TML estimate, accounting for censoring
tmle <- txshift(
  W = W, A = A, C_cens = C_cens, Y = Y, delta = 0.5,
  estimator = "onestep",
  g_exp_fit_args = list(
    fit_type = "hal",
    n_bins = 3,
    lambda_seq = exp(seq(-1, -10, length = 50))
  ),
  g_cens_fit_args = list(
    fit_type = "glm",
    glm_formula = "C_cens ~ ."
  ),
  Q_fit_args = list(
    fit_type = "glm",
    glm_formula = "Y ~ ."
  )
)

# construct a TML estimate under two-phase sampling, ignoring censoring
ipcwtmle <- txshift(
  W = W, A = A, Y = Y, delta = 0.5,
  C_samp = C_samp, V = c("W", "Y"),
  estimator = "onestep", max_iter = 3,
  samp_fit_args = list(fit_type = "glm"),
  g_exp_fit_args = list(
    fit_type = "hal",
    n_bins = 3,
    lambda_seq = exp(seq(-1, -10, length = 50))
  ),
  Q_fit_args = list(
    fit_type = "glm",
    glm_formula = "Y ~ ."
  ),
  eif_reg_type = "glm"
)

# construct a TML estimate acconting for two-phase sampling and censoring
ipcwtmle <- txshift(
  W = W, A = A, C_cens = C_cens, Y = Y, delta = 0.5,
  C_samp = C_samp, V = c("W", "Y"),
  estimator = "onestep", max_iter = 3,
  samp_fit_args = list(fit_type = "glm"),
  g_exp_fit_args = list(
    fit_type = "hal",
    n_bins = 3,
    lambda_seq = exp(seq(-1, -10, length = 50))
  ),
  g_cens_fit_args = list(
    fit_type = "glm",
    glm_formula = "C_cens ~ ."
  ),
  Q_fit_args = list(
    fit_type = "glm",
    glm_formula = "Y ~ ."
  ),
  eif_reg_type = "glm"
)
} # }
```
