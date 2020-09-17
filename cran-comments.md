## Test environments
* local Ubuntu 20.04: R 4.0.2 (stable)
* remote Ubuntu 16.04 on travis-ci:
  * R 3.6.3 (old-stable)
  * R 4.0.2 (stable)
* Windows (on appveyor and winbuilder): R 4.0.2

## R CMD check results
* There were 0 ERRORs.
* There were 0 WARNINGs.
* There were 0 NOTEs.

## Downstream dependencies
* Nothing to report.

## Additional Notes
* This is a new CRAN submission.
* In the initial attempt, there were two significant issues:
  * A few more links throughout used `http`; these have been moved to `https`
    or otherwise changed entirely.
  * The old URL `https://sl3.tlverse.org` has been replaced by
    `https://tlverse.org/sl3`
