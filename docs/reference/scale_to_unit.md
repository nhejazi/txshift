# Transform values by scaling to the unit interval

Transform values by scaling to the unit interval

## Usage

``` r
scale_to_unit(vals)
```

## Arguments

- vals:

  A `numeric` vector corresponding to the observed values of the
  variable of interest, to be re-scaled to the unit interval \[0,1\].

## Value

A `numeric` vector of the same length as `vals`, where the values are
re-scaled to lie in unit interval \[0, 1\].

## Details

A transformation that scales an arbitrary set of input values to the
unit interval. See
[`scale_to_original`](https://codex.nimahejazi.org/txshift/reference/scale_to_original.md)
for a corresponding backtransformation.
