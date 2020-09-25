# txshift 0.3.5

As of October 2020:
* Changes references to the argument `C` to `Delta` as the indicator of
  inclusion in the second-stage sample.
* Co-opts the newly freed argument `C` to denote censoring prior to occurrence
  of the outcome.
* Adds a nuisance regression for censoring `C` and adjusts the estimation
  procedure so as to use inverse censoring weights for this in the full-data
  EIF procedure (not the augmented two-phase sampling correction).

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
