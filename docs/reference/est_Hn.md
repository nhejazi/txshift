# Estimate Auxiliary Covariate of Full Data Efficient Influence Function

Estimate Auxiliary Covariate of Full Data Efficient Influence Function

## Usage

``` r
est_Hn(gn_exp, gps_bound)
```

## Arguments

- gn_exp:

  An estimate of the exposure density (a generalized propensity score)
  using the output provided by
  [`est_g_exp`](https://codex.nimahejazi.org/txshift/reference/est_g_exp.md).

- gps_bound:

  `numeric` giving the lower limit of the generalized propensity score
  estimates to be tolerated (default = 0.05). Estimates less than this
  are truncated to this or 1/n. See
  [`bound_propensity`](https://codex.nimahejazi.org/txshift/reference/bound_propensity.md)
  for details.

## Value

A `data.table` with two columns, containing estimates of the auxiliary
covariate at the natural value of the exposure H(A, W) and at the
shifted value of the exposure H(A + delta, W).

## Details

Compute an estimate of the auxiliary covariate of the efficient
influence function required to update initial estimates through logistic
tilting models for targeted minimum loss estimation.
