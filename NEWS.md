# txshift 0.3.9

As of August 2023:
* Txshift has been modified to estimate Q on the original scale of the outcome,
instead of scaling Y before estimating Q. This is helpful when the outcome is a 
count and some of the algorithms in the SL library only handle non-negative integers.
* A new argument `glm_family` has been added to `est_Q()` to specify the family
when a glm is used to estimate Q.
* @Lauren-EylerDang, #70

As of May 2023:
* A new argument `bound` has been added to `bound_propensity()` to specify the
  lower tolerated limit of generalized propensity score estimates. Estimates
  are bounded to the higher of the specified or default value of `bound` and
  the inverse of the sample size, 1/n.
* A new argument `gps_bound` has also been introduced to both `txshift()` and
  `est_Hn()` in order to accommodate passing in truncation bounds for the
  estimated generalized propensity score.

# txshift 0.3.6

As of October 2021:
* Minor updates to ensure compatibility with v0.4.1 of `hal9001` and v0.2.1 of
  `haldensify`, both recently updated on CRAN.
* Removal of the `LazyData` field from the `DESCRIPTION`, since no `data`
  directory is included with the package.
* Minor tweaks to existing unit tests to remove `rlang` from the `Suggests`
  field of the `DESCRIPTION`.
* Vignettes for the standard and IPCW-augmented estimation procedures have been
  combined to reduce redundancy and reduce build time per CRAN requests.

As of May 2021:
* The use of `hal9001::fit_hal()` internally for evaluation of a conditional
  mean of the full-data EIF has been revised for compatibility with v0.4.0+ of
  the `hal9001` package.
* Defaults passed in through the argument `g_exp_fit_args`, and to the function
  `est_g_exp()`, have been updated for compatibility with v0.1.5+ of the
  `haldensify` package.

As of April 2021:
* The `print()` methods have been updated to remove the use of [`cli`
  functions](https://github.com/r-lib/cli), which, for simplicity, has been
  replaced by the use of `message()`.
* Addition of a hidden slot `.eif_mat` to the `txshift_msm` class, supporting
  export of the matrix of EIF estimates for each shift in `delta_grid`.

# txshift 0.3.5

As of February 2021:
* Remove cross-linking to `sl3` functions as per request from CRAN. This can be
  reversed once `sl3` is available on CRAN.

As of January 2021:
* Simulation experiments testing the performance of the procedures in the
  presence of loss to follow-up censoring indicate that the TML estimator
  outperforms the one-step for the EIF-based two-phase sampling correction.
  Generally, we recommend use of the TML estimator (the default) across all
  settings, though performance of the one-step estimator is much worse.

As of December 2020:
* A `delta` slot has been added to the `txshift` class to record the shift.
* Hidden slots have been similarly added to the `txshift_msm` class.
* The `summary` method has been removed, with the functionality now supported
  by the `print` methods for the `txshift` and `txshift_msm` classes.
* The `plot` method has been amended to support simultaneous confidence bands.

As of October 2020:
* Changes all references to the argument `C` to `C_samp` for the indicator of
  inclusion in the second-stage sample.
* Adds the new argument `C_cens` to denote censoring due to loss to follow-up,
  i.e., prior to the occurrence of the outcome.
* Adds a nuisance regression for censoring `C_cens` and adjusts the estimation
  procedure so as to use inverse censoring weights in the full-data EIF
  procedure (NOTE: these are not updated in the two-phase sampling correction).
* Renaming of arguments to internal functions and functions themselves:
  * From `est_g` to `est_g_exp` for the exposure mechanism density estimation
  * From `est_ipcw` to `est_samp` for the two-phase sampling mechanism
  * Add `est_g_cens` for the loss to follow-up censoring mechanism

# txshift 0.3.4

As of September 2020:
* Moved `sl3` dependency to an `Enhances` designation for CRAN submission.
* As above, removed `sl3` from `Remotes` and added installation safety checks.

As of June 2020:
* Add single-knot spline to MSM summarization (`msm_vimshift`).
* Add class and `plot` method for MSM summarization (`msm_vimshift`).
* Fix bug in `msm_vimshift` for computing CIs for binary outcomes by switching
  from manually computing CIs to internally using custom `confint` method.
* Fix bug in `msm_vimshift` for building `lm` model objects through weighted
  regression; move models from `plot` method to `msm_vimshift`.
* Finish drafting brief paper for _Journal of Open Source Software_.

# txshift 0.3.3

As of April 2020:
* Change export status of internal functions (e.g., no longer exporting
  `onestep_txshift` and `tmle_txshift`).
* Finish adding Roxygen "details" and "return" slots throughout functions.
* Add examples to main estimation functions (`txshift`, `vimshift_msm`).
* Update argument names and add several `assert_that` checks.
* Change `fit_spec` terminology to `fit_ext` for external fits.
* Add unit tests for MSM functionality and nuisance parameter estimation.

As of March 2020:
* Extensive documentation, including fixing estimation terminology (e.g.,
  one-step instead of AIPW) and adding Roxygen "details" and "return" slots.
* Begin adding examples to exported functions.

# txshift 0.3.2

As of March 2020:
* Corrections to dependencies in preparation for eventual CRAN release.
* Change several previously exported functions to internal, including `eif`,
  `est_Hn`, `est_Q`, `est_g`, `est_ipcw`, `fit_fluctuation`, `ipcw_eif_update`).
* Remove/reduce GitHub-only dependencies (now only `sl3`).
* Change title partially (from "Targeted Learning" to "Efficient Estimation").
* Lock dependency versions (e.g., `sl3` >= v1.3.7)
* Extensive documentation updates.

# txshift 0.3.1

As of December 2019:
* Changes arguments of `hal9001::fit_hal` in pseudo-outcome regression for
    efficient estimation by explicitly including `max_degree = NULL`.
* Change to TMLE convergence criterion: use a less strict criterion such that
     | Pn D | \leq sigma / (sqrt(n) \cdot max(10, log(n))) instead of \leq 1/n.
    Empirical studies suggest this curbs issues addressed by over-agressive
    updates from the targeting step.
* Remove pinning of `sl3` dependency to a specific tag (formerly v1.2.0).
* Lock dependency version: `sl3` >= v1.3.6 and `hal9001` >= v0.2.5.

# txshift 0.3.0

As of October 2019:
* Change use of `as.data.table` to `data.table` in internal functions to catch
    up with changes in dependencies.

# txshift 0.2.9

As of September 2019:
* Remove errant intercept term and lower iterations for fluctuation models.
* Change weighting scheme in marginal structural model summarization to weight
    all estimates identically rather than by inverse variance as a default.
* Updates to documentation.

# txshift 0.2.8

As of September 2019:
* Add safety checks for convergence of fluctuation regressions based on those
    appearing in `drtmle` and/or `survtmle`.
* Change default confidence interval type to use marginal CIs across multiple
    parameters instead of a simultaneous confidence band.
* Switch internal parametric regressions to use `sl3::Lrnr_glm` instead of
    `sl3::Lrnr_glm_fast`.

# txshift 0.2.7

As of July 2019:
* Improve argument names for clarity and update documentation.
* Addition of tighter unit tests for both one-step and TML estimators.

As of June 2019:
* Pin `sl3` dependency to version 1.2.0 of that package for stability.

# txshift 0.2.6

As of June 2019:
* Changes to arguments of `hal9001::fit_hal` for pseudo-outcome EIF regression.
* Addition of clarifying notes to core internal functions.
* Removal of outdated (and commented out) code in core internal functions.
* Clarifying alterations to internal function and argument names.
* Renaming internal function `tx_shift` to `shift_additive`.

# txshift 0.2.5

As of June 2019:
* Remove inverse weights from estimated efficient influence function necessary
    for pseudo-outcome regression for efficient IPCW-augmented estimators.
* Reduce use of redundant variables across core functions, reorganize functions
    across files, clarify documentation.
* Tweak arguments for fitting pseudo-outcome regression with HAL in order to
    diagnose performance issues revealed by simulation.
* Fix how inverse weights are passed to full-data estimators.
* Pare down arguments for the one-step estimation routine.

# txshift 0.2.4

As of June 2019:
* Introduce `bound_propensity` function to bound the propensity score away
    from zero by a factor 1/n, rather than to numerical precision.

As of April 2019:
* Minor improvements to documentation and vignettes.
* Fix a bug in the output of the IPCW one-step estimator.
* Pare down packages listed in imports, moving several to suggests.
* Introduce option to compute simultaneous confidence band for working MSMs.
* Fix a bug introduced by newly added imputation functionality in `sl3`.

# txshift 0.2.3

As of March 2019:
* Introduce functionality for computing one-step estimators to complement the
    the available TMLEs.
* Add initial functionality for summarizing estimated effects across a grid of
    shifts via working marginal structural models (MSMs).

# txshift 0.2.2

As of February 2019:
* Added helper functions and caught edge cases in auxiliary covariate for TMLE
    fluctuation models.
* Fixed a bug in how the auxiliary covariate for TMLEs is computed by keeping
    track of an extra shift g(a+2*delta|w).
* Revised inference machinery to create confidence intervals on the logit scale
    in the case of binary outcomes.

# txshift 0.2.0

As of May 2018:
* An initial public release of this package, version 0.2.0.
* This version including complete functionality for both standard TML and
    IPCW-TML estimators.
