# Estimate the Exposure Mechanism via Generalized Propensity Score

Estimate the Exposure Mechanism via Generalized Propensity Score

## Usage

``` r
est_g_exp(
  A,
  W,
  delta = 0,
  samp_weights = rep(1, length(A)),
  fit_type = c("hal", "sl"),
  sl_learners_density = NULL,
  haldensify_args = list(grid_type = "equal_range", lambda_seq = exp(seq(-1, -13, length
    = 300)))
)
```

## Arguments

- A:

  A `numeric` vector of observed exposure values.

- W:

  A `numeric` matrix of observed baseline covariate values.

- delta:

  A `numeric` value identifying a shift in the observed value of the
  exposure under which observations are to be evaluated.

- samp_weights:

  A `numeric` vector of observation-level sampling weights, as produced
  by the internal procedure to estimate the two-phase sampling mechanism
  [`est_samp`](https://codex.nimahejazi.org/txshift/reference/est_samp.md).

- fit_type:

  A `character` specifying whether to use Super Learner (from sl3) or
  the Highly Adaptive Lasso (from hal9001) to estimate the conditional
  exposure density.

- sl_learners_density:

  Object containing a set of instantiated learners from sl3, to be used
  in fitting an ensemble model.

- haldensify_args:

  A `list` of argument to be directly passed to
  [`haldensify`](https://codex.nimahejazi.org/haldensify/reference/haldensify.html)
  when `fit_type` is set to `"hal"`. Note that this invokes the Highly
  Adaptive Lasso instead of Super Learner and is thus only feasible for
  relatively small data sets.

## Value

A `data.table` with four columns, containing estimates of the
generalized propensity score at a downshift (g(A - delta \| W)), no
shift (g(A \| W)), an upshift (g(A + delta) \| W), and an upshift of
magnitude two (g(A + 2 delta) \| W).

## Details

Compute the propensity score (exposure mechanism) for the observed data,
including the shift. This gives the propensity score for the observed
data (at the observed A) the counterfactual shifted exposure levels at
(A - delta), (A + delta), and (A + 2 \* delta).
