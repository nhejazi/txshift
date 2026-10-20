# Estimate Probability of Censoring by Two-Phase Sampling

Estimate Probability of Censoring by Two-Phase Sampling

## Usage

``` r
est_samp(V, C_samp, fit_type = c("sl", "glm"), sl_learners = NULL)
```

## Arguments

- V:

  A `numeric` vector, `matrix`, `data.frame` or similar object giving
  the observed values of the covariates known to potentially inform the
  sampling mechanism.

- C_samp:

  A `numeric` vector of observed values of the indicator for inclusion
  in the second-phase sample.

- fit_type:

  A `character` indicating whether to perform the fit using GLMs or a
  Super Learner ensemble model. If use of Super Learner is desired, then
  the argument `sl_learners` must be provided.

- sl_learners:

  An sl3 `Lrnr_sl` object, a Super Learner ensemble or learner
  instantiated externally using sl3.

## Value

A `numeric` vector of the estimated sampling mechanism.

## Details

Compute estimates of the sampling probability for inclusion in the the
second-phase via the two-phase sampling mechanism. These estimates are
used for the creation of inverse probability weights.
