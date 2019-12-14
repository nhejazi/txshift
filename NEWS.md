# txshift 0.3.1

As of December 2019:
* ...
* ...

# txshift 0.3.0

As of October 2019:
* ...
* ...

# txshift 0.2.9

As of September 2019:
* ...
* ...

# txshift 0.2.8

As of September 2019:
* ...
* ...

# txshift 0.2.7

As of July 2019:
* ...
* ...

As of June 2019:
* ...
* ...

# txshift 0.2.6

As of June 2019:
* ...
* ...

# txshift 0.2.5

As of June 2019:
* ...
* ...

# txshift 0.2.4

As of June 2019:
* ...
* ...

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
