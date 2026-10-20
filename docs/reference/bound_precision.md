# Bound Precision

Bound Precision

## Usage

``` r
bound_precision(vals)
```

## Arguments

- vals:

  `numeric` vector of values in the interval \[0, 1\] to be bounded
  within arbitrary machine precision. The most common use of this
  functionality is to avoid indeterminate or non-finite values after the
  application [`stats::qlogis`](https://rdrr.io/r/stats/Logistic.html).

## Value

A `numeric` vector of the same length as `vals`, where the returned
values are bounded to machine precision. This is intended to avoid
numerical instability issues.

## Details

Bound values in the unit interval to machine precision in order to avoid
numerical instability issues in downstream computation.
