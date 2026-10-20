# Confidence Intervals for Counterfactual Mean Under Stochastic Intervention

Confidence Intervals for Counterfactual Mean Under Stochastic
Intervention

## Usage

``` r
# S3 method for class 'txshift'
confint(object, parm = seq_len(object$psi), level = 0.95, ..., ci_mult = NULL)
```

## Arguments

- object:

  An object of class `txshift`, as produced by invoking the function
  [`txshift`](https://codex.nimahejazi.org/txshift/reference/txshift.md),
  for which a confidence interval is to be computed.

- parm:

  A `numeric` vector indicating indices of `object$est` for which to
  return confidence intervals.

- level:

  A `numeric` indicating the level of the confidence interval to be
  computed.

- ...:

  Other arguments. Not currently used.

- ci_mult:

  Pre-computed multipliers for generating confidence intervals. The
  default of `NULL` should generally NOT be changed and is only used by
  the internal machinery for creating simultaneous confidence bands.

## Value

A named `numeric` vector containing the parameter estimate from a
`txshift` object, alongside lower and upper Wald-style confidence
intervals at a specified coverage level.

## Details

Compute confidence intervals for estimates produced by
[`txshift`](https://codex.nimahejazi.org/txshift/reference/txshift.md).

## Examples

``` r
set.seed(429153)
n_obs <- 100
W <- replicate(2, rbinom(n_obs, 1, 0.5))
A <- rnorm(n_obs, mean = 2 * W, sd = 1)
Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
txout <- txshift(
  W = W, A = A, Y = Y, delta = 0.5,
  estimator = "tmle",
  g_exp_fit_args = list(
    fit_type = "hal", n_bins = 5,
    grid_type = "equal_mass",
    lambda_seq = exp(-1:-9)
  ),
  Q_fit_args = list(
    fit_type = "glm",
    glm_formula = "Y ~ ."
  )
)
confint(txout)
#>    lwr_ci       est    upr_ci 
#> 0.6918278 0.7845510 0.8552141 
```
