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
* This is a re-submission of a CRAN update that fixed an issue raised by CRAN:
  * Instances of `\code{\link[sl3]{...}}` in the documentation were changed to
      remove the `\link[]` statement referencing a package unavailable on CRAN.
  * This re-submission fixes an issue concerning http vs. https links:
                    Found the following (possibly) invalid URLs:
                    URL: http://www.r-pkg.org/pkg/txshift (moved to
                      https://www.r-pkg.org/pkg/txshift)
                             From: README.md
                             Status: 200
                             Message: OK
