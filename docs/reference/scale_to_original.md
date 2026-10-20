# Transform values from the unit interval back to their original scale

Transform values from the unit interval back to their original scale

## Usage

``` r
scale_to_original(scaled_vals, max_orig, min_orig)
```

## Arguments

- scaled_vals:

  A `numeric` vector corresponding to re-scaled values in the unit
  interval, to be re-scaled to the original interval.

- max_orig:

  A `numeric` scalar value giving the maximum of the values on the
  original scale.

- min_orig:

  A `numeric` scalar value giving the minimum of the values on the
  original scale.

## Value

A `numeric` vector of the same length as `scaled_vals`, where the values
are re-scaled to lie in their original/natural interval.

## Details

A back-transformation that returns values computed in the unit interval
to their original scale. This is used in re-scaling updated TML
estimates back to their natural scale. Undoes
[`scale_to_unit`](https://codex.nimahejazi.org/txshift/reference/scale_to_unit.md).
