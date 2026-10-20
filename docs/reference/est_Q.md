# Estimate the Outcome Mechanism

Estimate the Outcome Mechanism

## Usage

``` r
est_Q(
  Y,
  C_cens = rep(1, length(Y)),
  A,
  W,
  delta = 0,
  samp_weights = rep(1, length(Y)),
  fit_type = c("sl", "glm"),
  glm_formula = "Y ~ .",
  sl_learners = NULL,
  scale_outcome = TRUE
)
```

## Arguments

- Y:

  A `numeric` vector of observed outcomes.

- C_cens:

  A `numeric` vector of loss to follow-up indicators.

- A:

  A `numeric` vector of observed exposure values.

- W:

  A `numeric` matrix of observed baseline covariate values.

- delta:

  A `numeric` indicating the magnitude of the shift to be computed for
  the exposure `A`. This is passed to the internal
  [`shift_additive`](https://codex.nimahejazi.org/txshift/reference/shift_additive.md)
  and is currently limited to additive shifts.

- samp_weights:

  A `numeric` vector of observation-level sampling weights, as produced
  by the internal procedure to estimate the two-phase sampling mechanism
  [`est_samp`](https://codex.nimahejazi.org/txshift/reference/est_samp.md).

- fit_type:

  A `character` indicating whether to use GLMs or Super Learner to fit
  the outcome regression. If the option "glm" is selected, the argument
  `glm_formula` must NOT be `NULL`, instead containing a model formula
  (as per [`glm`](https://rdrr.io/r/stats/glm.html)) as a `character`.
  If the option "sl" is selected, the argument `sl_learners` must NOT be
  `NULL`; instead, an instantiated sl3 `Lrnr_sl` object, specifying
  learners and a metalearner for the Super Learner fit, must be
  provided. Consult the documentation of sl3 for details.

- glm_formula:

  A `character` giving a
  [`formula`](https://rdrr.io/r/stats/formula.html) for fitting a
  (generalized) linear model via
  [`glm`](https://rdrr.io/r/stats/glm.html).

- sl_learners:

  Object containing a set of instantiated learners from the sl3, to be
  used in fitting an ensemble model.

- scale_outcome:

  A `logical` indicating whether the outcome should be scaled to the
  unit interval before estimating the outcome regression (the default).
  If `FALSE`, the outcome regression is estimated on the original scale
  of the outcome (e.g., for count outcomes, with learners that require
  non-negative integers), and its predictions are then scaled to the
  unit interval, as required by the targeting step. When `FALSE`,
  GLM-based fits use a Gaussian family; other working models (e.g.,
  Poisson) may be specified via `sl_learners`.

## Value

A `data.table` with two columns, containing estimates of the outcome
mechanism at the natural value of the exposure Q(A, W) and an upshift of
the exposure Q(A + delta, W).

## Details

Compute the outcome regression for the observed data, including with the
shift imposed by the intervention. This returns the outcome regression
for the observed data (at A) and under the counterfactual shift shift
(at A + delta).
