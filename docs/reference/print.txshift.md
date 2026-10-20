# Print Method for Counterfactual Mean of Stochastic Shift Intervention

Print Method for Counterfactual Mean of Stochastic Shift Intervention

## Usage

``` r
# S3 method for class 'txshift'
print(x, ..., ci_level = 0.95)
```

## Arguments

- x:

  An object of class `txshift`.

- ...:

  Other options (not currently used).

- ci_level:

  A `numeric` indicating the level of the confidence interval to be
  computed.

## Value

None. Called for the side effect of printing an informative summary of
slots of objects of class `txshift`.

## Details

The `print` method for objects of class `txshift`.

## Examples

``` r
if (require("sl3")) {
  set.seed(429153)
  n_obs <- 100
  W <- replicate(2, rbinom(n_obs, 1, 0.5))
  A <- rnorm(n_obs, mean = 2 * W, sd = 1)
  Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))
  txout <- txshift(
    W = W, A = A, Y = Y, delta = 0.5,
    estimator = "tmle",
    g_exp_fit_args = list(
      fit_type = "sl",
      sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
    ),
    Q_fit_args = list(
      fit_type = "glm",
      glm_formula = "Y ~ ."
    )
  )
  print(txout)
}
#> Counterfactual Mean of Shifted Treatment
#> Intervention: Treatment + 0.5
#> txshift Estimator: tmle
#> Estimate: 0.7876
#> Std. Error: 0.0387
#> 95% CI: [0.702, 0.8538]
```
