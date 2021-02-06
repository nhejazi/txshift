## Test environments
* local Ubuntu 20.04: R 4.0.3 (stable)
* remote Ubuntu 18.04 on Travis-CI:
  * R 3.6.3 (old-stable)
  * R 4.0.3 (stable)
* Windows (on Appveyor and Winbuilder): R 4.0.3

## R CMD check results
* There were 0 ERRORs.
* There were 0 WARNINGs.
* There were 0 NOTEs.

## Downstream dependencies
* Nothing to report.

## Additional Notes
* This is an updated CRAN submission to fix an issue raised by CRAN.
* The issue fixed concerns cross-referencing a package not available on CRAN:
  * Instances of `\code{\link[sl3]{...}}` in the documentation were changed to
      remove the `\link[]` statement.
