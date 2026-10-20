# One-Step Estimate of Counterfactual Mean of Stochastic Shift Intervention

One-Step Estimate of Counterfactual Mean of Stochastic Shift
Intervention

## Usage

``` r
onestep_txshift(
  data_internal,
  C_samp = rep(1, nrow(data_internal)),
  V = NULL,
  delta,
  samp_estim,
  gn_cens_weights,
  Qn_estim,
  Hn_estim,
  eif_reg_type = c("hal", "glm"),
  samp_fit_args,
  ipcw_efficiency = TRUE
)
```

## Arguments

- data_internal:

  A `data.table` constructed internally by a call to
  [`txshift`](https://codex.nimahejazi.org/txshift/reference/txshift.md).
  This contains the data elements needed for computing the one-step
  estimator.

- C_samp:

  A `numeric` indicator for whether a given observation was included in
  the second-stage sample, used to compute an IPC-weighted one-step
  estimator in cases where two-stage sampling is performed. Default
  assumes no censoring due to sampling.

- V:

  The covariates that are used in determining the sampling procedure
  that gives rise to censoring. The default is `NULL` and corresponds to
  scenarios in which there is no censoring (in which case all values in
  the preceding argument `C` must be uniquely 1. To specify this, pass
  in a NAMED `list` identifying variables amongst W, A, Y that are
  thought to have played a role in defining the sampling/censoring
  mechanism (C).

- delta:

  A `numeric` value indicating the shift in the treatment to be used in
  defining the target parameter. This is defined with respect to the
  scale of the treatment (A).

- samp_estim:

  An object providing the value of the censoring mechanism evaluated
  across the full data. This object is passed in after being constructed
  by a call to the internal function
  [`est_samp`](https://codex.nimahejazi.org/txshift/reference/est_samp.md).

- gn_cens_weights:

  An object providing the value of inverse probability of censoring
  weights, the inverse of the censoring mechanism estimate. The weights
  are used as part of the IPCW-EIF procedure to implement a joint
  intervention that removes the contribution of the censoring process.

- Qn_estim:

  An object providing the value of the outcome evaluated after imposing
  a shift in the treatment. This object is passed in after being
  constructed by a call to the internal function `est_Q`.

- Hn_estim:

  An object providing values of the auxiliary ("clever") covariate,
  constructed from the treatment mechanism and required for targeted
  minimum loss estimation. This object object should be passed in after
  being constructed by a call to the internal function `est_Hn`.

- eif_reg_type:

  Whether a flexible nonparametric function ought to be used in the
  dimension-reduced nuisance regression of the targeting step for the
  censored data case. By default, the method used is a nonparametric
  regression based on the Highly Adaptive Lasso (from hal9001). Set this
  to `"glm"` to instead use a simple linear regression model. In this
  step, the efficient influence function (EIF) is regressed against
  covariates contributing to the censoring mechanism (i.e., EIF ~ V \| C
  = 1).

- samp_fit_args:

  A `list` of arguments, all but one of which are passed to
  [`est_samp`](https://codex.nimahejazi.org/txshift/reference/est_samp.md).
  For details, consult the documentation for
  [`est_samp`](https://codex.nimahejazi.org/txshift/reference/est_samp.md).
  The first element (i.e., `fit_type`) is used to determine how this
  regression is fit: "glm" for generalized linear model, "sl" for a
  Super Learner, and "external" for a user-specified input of the form
  produced by
  [`est_samp`](https://codex.nimahejazi.org/txshift/reference/est_samp.md).

- ipcw_efficiency:

  Whether to invoke an augmentation of the IPCW-TMLE procedure that
  performs an iterative process to ensure efficiency of the resulting
  estimate. The default is `TRUE`; set to `FALSE` to use an IPC-weighted
  loss rather than the IPC-augmented influence function.

## Value

S3 object of class `txshift` containing the results of the procedure to
compute a one-step estimate of the treatment shift parameter.

## Details

Invokes the procedure to construct a one-step estimate of the
counterfactual mean under a modified treatment policy.
