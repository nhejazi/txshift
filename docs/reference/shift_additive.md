# Simple Additive Modified Treatment Policy

Simple Additive Modified Treatment Policy

## Usage

``` r
shift_additive(A, W = NULL, delta)
```

## Arguments

- A:

  A `numeric` vector of observed treatment values.

- W:

  A `numeric` matrix of observed baseline covariate values.

- delta:

  A `numeric` indicating the magnitude of the shift to be computed for
  the treatment `A`.

## Value

A `numeric` vector containing the shifted exposure values.

## Details

A simple modified treatment policy that modifes the observed value of
the exposure by shifting it by a value `delta`. Note that this shifting
function assumes support of A\|W across all strata of W.
