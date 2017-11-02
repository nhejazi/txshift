
<!-- README.md is generated from README.Rmd. Please edit that file -->
R/`shifttx`
===========

[![Travis-CI Build
Status](https://travis-ci.org/nhejazi/shifttx.svg?branch=master)](https://travis-ci.org/nhejazi/shifttx)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/nhejazi/shifttx?branch=master&svg=true)](https://ci.appveyor.com/project/nhejazi/shifttx)
[![Coverage
Status](https://img.shields.io/codecov/c/github/nhejazi/shifttx/master.svg)](https://codecov.io/github/nhejazi/shifttx?branch=master)
[![Project Status: WIP - Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> …

**Authors:** [Nima Hejazi](http://nimahejazi.org) and [David
Benkeser](https://www.benkeserstatistics.com/)

------------------------------------------------------------------------

What’s `shifttx`?
-----------------

`shifttx` is an R package that makes it easy to compute targeted minimum
loss-based estimates of the population-level causal effects of
interventions based on stochastic mechanisms for treatment assignment.

-   Original estimator and iterative algorithm proposed in Muñoz and van
    der Laan (2012)
-   One-step estimation procedure and algorithm introduced in Díaz and
    van der Laan (2017).
-   See van der Laan and Rose (2017) for a discussion on recent
    developments in Targeted Learning.
-   See van der Laan and Rose (2011) for an introduction to Targeted
    Learning.

------------------------------------------------------------------------

Installation
------------

Install the most recent *stable release* from GitHub via
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/):

``` r
devtools::install_github("nhejazi/shifttx")
```

------------------------------------------------------------------------

Example
-------

To illustrate how `shifttx` may be used to ascertain the effect of a
treatment, consider the following example:

``` r
library(shifttx)
#> shifttx: Estimate Causal Effects with Stochastic Treatments
#> Version: 0.0.0.9000

# first, define simulation parameters and observed data
n <- 100
W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.7))
A <- rpois(n, lambda = exp(3 + 0.3 * log(W$W1) - 0.2 * exp(W$W1) * W$W2))
Y <- rbinom(n, 1, plogis(-1 + 0.05 * A - 0.02 * A * W$W2 + 0.2 * A * tan(W$W1^2)
                         - 0.02 * W$W1 * W$W2 + 0.1 * A * W$W1 * W$W2))
fitA.0 <- glm(A ~ I(log(W1)) + I(exp(W1)):W2, family = poisson,
              data = data.frame(A, W))
fitY.0 <- glm(Y ~ A + A:W2 + A:I(tan(W1^2)) + W1:W2 + A:W1:W2,
              family = binomial, data = data.frame(A, W))

# next, we define the treatment mechanisma and outcome in the true model
gn.0  <- function(A = A, W = W) {
  dpois(A, lambda = predict(fitA.0, newdata = W, type = "response"))
}
Qn.0 <- function(A = A, W = W) {
  predict(fitY.0, newdata = data.frame(A, W, row.names = NULL),
          type = "response")
}

# finally, we may compute targeted estimates of the intervention effect
est_target <- tmle_shift(Y = Y, A = A, W = W, Qn = Qn.0, gn = gn.0, delta = 2,
                         tol = 1e-4, iter_max = 5, A_val = seq(1, 60, 1))

# let's just examine the first few results
head(est_target)
#>      psi.hat      var.hat         IC.1         IC.2         IC.3 
#>  0.584564153  0.003225504  0.429013014 -1.007698967  0.462008069 
#>         IC.4 
#>  0.420323583
```

------------------------------------------------------------------------

Issues
------

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/nhejazi/shifttx/issues).

------------------------------------------------------------------------

License
-------

© 2017 [Nima S. Hejazi](http://nimahejazi.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License

    Copyright (c) 2017 Nima S. Hejazi

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

------------------------------------------------------------------------

References
----------

Díaz, Iván, and Mark J van der Laan. 2017. “Stochastic Treatment
Regimes.” In *Targeted Learning in Data Science: Causal Inference for
Complex Longitudinal Studies*, 167–80. Springer Science & Business
Media.

Muñoz, Iván Díaz, and Mark J van der Laan. 2012. “Population
Intervention Causal Effects Based on Stochastic Interventions.”
*Biometrics* 68 (2). Wiley Online Library: 541–49.

van der Laan, Mark J, and Sherri Rose. 2011. *Targeted Learning: Causal
Inference for Observational and Experimental Data*. Springer Science &
Business Media.

———. 2017. *Targeted Learning in Data Science: Causal Inference for
Complex Longitudinal Studies*. Springer Science & Business Media.
