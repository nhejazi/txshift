# Bound Generalized Propensity Score

Bound Generalized Propensity Score

## Usage

``` r
bound_propensity(vals, bound = 0.01)
```

## Arguments

- vals:

  `numeric` vector of generalized propensity score estimates. Note that
  the generalized propensity score is a (conditional) density and should
  only be bounded/truncated away from zero.

- bound:

  `numeric` atomic giving the lower limit of the generalized propensity
  score estimates to be tolerated. Estimates less than this will be
  truncated to this value (default = 0.01). Note that the default may be
  internally overwritten by reference to the sample size (n), so the
  lower bound used is the greater of the specified value and 1/n.

## Value

A `numeric` vector of the same length as `vals`, where the returned
values are bounded such that the minimum is no lower than 1/n, for the
sample size n.

## Details

Bound estimated values of the generalized propensity score (a
conditional density) to avoid numerical instability issues arising from
practical violations of the assumption of positivity.
