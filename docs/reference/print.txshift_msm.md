# Print Method for Marginal Structural Models

Print Method for Marginal Structural Models

## Usage

``` r
# S3 method for class 'txshift_msm'
print(x, ...)
```

## Arguments

- x:

  An object of class `txshift_msm`.

- ...:

  Other options (not currently used).

## Value

None. Called for the side effect of printing an informative summary of
slots of objects of class `txshift_msm`.

## Details

The `print` method for objects of class `txshift_msm`.

## Examples

``` r
if (require("sl3")) {
  set.seed(3287)
  n_obs <- 1000
  W <- as.numeric(replicate(1, rbinom(n_obs, 1, 0.5)))
  A <- as.numeric(rnorm(n_obs, mean = 2 * W, sd = 1))
  Y <- rbinom(n_obs, 1, plogis(2 * A - W))
  msm <- msm_vimshift(
    W = W, A = A, Y = Y, estimator = "tmle",
    g_exp_fit_args = list(
      fit_type = "sl",
      sl_learners_density = Lrnr_density_hse$new(Lrnr_glm$new())
    ),
    Q_fit_args = list(
      fit_type = "glm",
      glm_formula = "Y ~ ."
    ),
    delta_grid = seq(-1, 1, 0.25)
  )
  print(msm)
}
#> MSM (linear) for Grid of Shifted Treatments
#> Intervention Grid: Treatment + {-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1}
#> txshift MSM Estimator: tmle
#> Estimated Slope: 0.2114
#> Std. Error: 0.0095
#> 95% CI: [0.1928, 0.2301]
#> p-value (vs. no trend): 0
```
