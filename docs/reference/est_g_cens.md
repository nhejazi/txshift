# Estimate the Censoring Mechanism

Estimate the Censoring Mechanism

## Usage

``` r
est_g_cens(
  C_cens,
  A,
  W,
  samp_weights = rep(1, length(C_cens)),
  fit_type = c("sl", "glm"),
  glm_formula = "C_cens ~ .",
  sl_learners = NULL,
  bound = 0.02
)
```

## Arguments

- C_cens:

  A `numeric` vector of loss to follow-up indicators.

- A:

  A `numeric` vector of observed exposure values.

- W:

  A `numeric` matrix of observed baseline covariate values.

- samp_weights:

  A `numeric` vector of observation-level sampling weights, as produced
  by the internal procedure to estimate the two-phase sampling mechanism
  [`est_samp`](https://codex.nimahejazi.org/txshift/reference/est_samp.md).

- fit_type:

  A `character` indicating whether to use GLMs or Super Learner to fit
  the censoring mechanism. If option "glm" is selected, the argument
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

- bound:

  `numeric` giving the lower limit of censoring mechanism estimates to
  be tolerated (default = 0.02). Estimates below this value are
  truncated to this or 1/n. See
  [`bound_propensity`](https://codex.nimahejazi.org/txshift/reference/bound_propensity.md)
  for details.

## Value

A `numeric` vector of the propensity score for censoring.

## Details

Compute the censoring mechanism for the observed data, in order to apply
a joint intervention for removing censoring by re-weighting.
