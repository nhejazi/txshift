# Iterative IPCW Update Procedure of Augmented Efficient Influence Function

Iterative IPCW Update Procedure of Augmented Efficient Influence
Function

## Usage

``` r
ipcw_eif_update(
  data_internal,
  C_samp,
  V,
  ipc_mech,
  ipc_weights,
  Qn_estim,
  Hn_estim,
  estimator = c("tmle", "onestep"),
  fluctuation = NULL,
  flucmod_tol = 50,
  eif_reg_type = c("hal", "glm")
)
```

## Arguments

- data_internal:

  A `data.table` containing of the observations selected into the
  second-phase sample.

- C_samp:

  A `numeric` indicator for missingness due to exclusion the from
  second-stage sample.

- V:

  A `data.table` giving the values across all observations of all
  variables that play a role in the censoring mechanism.

- ipc_mech:

  A `numeric` vector of the censoring mechanism estimates all of the
  observations, only for the two-phase sampling mechanism. Note well
  that these values do NOT account for censoring from loss to follow-up.

- ipc_weights:

  A `numeric` vector of inverse probability of censoring weights,
  including such weights for censoring due to loss to follow-up. Without
  loss to follow-up, these are equivalent to `C_samp / ipc_mech` in an
  initial run of this procedure.

- Qn_estim:

  A `data.table` corresponding to the outcome regression. This is
  produced by invoking the internal function `est_Q`.

- Hn_estim:

  A `data.table` corresponding to values produced in the computation of
  the auxiliary ("clever") covariate. This is produced easily by
  invoking the internal function `est_Hn`.

- estimator:

  The type of estimator to be fit, either `"tmle"` for targeted maximum
  likelihood estimation or `"onestep"` for a one-step estimator.

- fluctuation:

  A `character` giving the type of regression to be used in traversing
  the fluctuation submodel. The choices are "weighted" and "standard".

- flucmod_tol:

  A `numeric` indicating the largest value to be tolerated in the
  fluctuation model for the targeted minimum loss estimator.

- eif_reg_type:

  Whether a flexible nonparametric function ought to be used in the
  dimension-reduced nuisance regression of the targeting step for the
  censored data case. By default, the method used is a nonparametric
  regression based on the Highly Adaptive Lasso (from hal9001). Set this
  to `"glm"` to instead use a simple linear regression model. In this
  step, the efficient influence function (EIF) is regressed against
  covariates contributing to the censoring mechanism (i.e., EIF ~ V \| C
  = 1).

## Value

A `list` containing the estimated outcome mechanism, the fitted
fluctuation model for TML updates, the updated inverse probability of
censoring weights (IPCW), the updated estimate of the efficient
influence function, and the estimated IPCW component of the EIF.

## Details

An adaptation of the IPCW-TMLE for iteratively constructing an efficient
inverse probability of censoring weighted TML or one-step estimator. The
efficient influence function of the parameter and updating the IPC
weights in an iterative process, until a convergence criteria is
satisfied.
