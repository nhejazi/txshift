
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`txshift`

[![Travis-CI Build
Status](https://travis-ci.org/nhejazi/txshift.svg?branch=master)](https://travis-ci.org/nhejazi/txshift)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/nhejazi/txshift?branch=master&svg=true)](https://ci.appveyor.com/project/nhejazi/txshift)
[![Coverage
Status](https://img.shields.io/codecov/c/github/nhejazi/txshift/master.svg)](https://codecov.io/github/nhejazi/txshift?branch=master)
[![Project Status: WIP - Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> Targeted Learning of Continuous Intervention Effects with Stochastic
> Treatment Regimes

**Authors:** [Nima Hejazi](http://nimahejazi.org) and [David
Benkeser](https://www.benkeserstatistics.com/)

-----

## What’s `txshift`?

`txshift` is an R package that makes it easy to compute targeted minimum
loss-based estimates of the population-level causal effects of
interventions based on stochastic mechanisms for treatment assignment.

  - Original estimator and iterative algorithm proposed in Muñoz and van
    der Laan (2012)
  - One-step estimation procedure and algorithm introduced in Díaz and
    van der Laan (2017).
  - See van der Laan and Rose (2017) for a discussion on recent
    developments in Targeted Learning.
  - See van der Laan and Rose (2011) for an introduction to Targeted
    Learning.

-----

## Installation

Install the most recent *stable release* from GitHub via
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/):

``` r
devtools::install_github("nhejazi/txshift", build_vignettes = FALSE)
```

-----

## Example

To illustrate how `txshift` may be used to ascertain the effect of a
treatment, consider the following example:

``` r
library(tidyverse)
#> ── Attaching packages ────────────────────────────────── tidyverse 1.2.1 ──
#> ✔ ggplot2 2.2.1.9000     ✔ purrr   0.2.4     
#> ✔ tibble  1.4.2          ✔ dplyr   0.7.4     
#> ✔ tidyr   0.7.2.9000     ✔ stringr 1.2.0.9000
#> ✔ readr   1.1.1          ✔ forcats 0.2.0
#> ── Conflicts ───────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
library(condensier)
#> condensier
#> The condensier package is still in beta testing. Interpret results with caution.
library(txshift)
#> txshift: Targeted Learning with Stochastic Interventions
#> Version: 0.0.9.3
set.seed(429153)

# simulate simple data for tmle-shift sketch
n_obs <- 1000  # number of observations
n_w <- 1  # number of baseline covariates
p_w <- 0.5  # probability of a success ("1") in the baseline variables
tx_mult <- 2  # multiplier for the effect of W = 1 on the treatment

## baseline covariate -- simple, binary
W <- as.numeric(replicate(n_w, rbinom(n_obs, 1, p_w)))

## create treatment based on baseline W
A <- as.numeric(rnorm(n_obs, mean = tx_mult * W, sd = 1))

# create outcome as a linear function of A, W + white noise
Y <- A + W + rnorm(n_obs, mean = 0, sd = 1)

# induce censoring
C <- rbinom(n_obs, 1, plogis(W))

# fit the TMLE
tmle_shift <- tmle_txshift(W = W, A = A, Y = Y, delta = 0.5,
                           fluc_method = "standard",
                           mod_args = list(
                             g_fit = list(fit_type = "glm",
                                          nbins = 25,  bin_method = "dhist",
                                          bin_estimator = speedglmR6$new(),
                                          parfit = FALSE),
                             Q_fit = list(fit_type = "glm",
                                          glm_formula = "Y ~ .")
                          ))

# examine the results
tmle_shift
#> $psi
#> [1] 2.060616
#> 
#> $var
#> [1] 0.004902471
#> 
#> $msg
#> [1] "EIF mean < 1e-09 (sufficiently low)."

# compute the confidence interval and view the results
(ci_shift <- confint(tmle_shift))
#>   lwr_CI      est   upr_CI 
#> 1.923384 2.060616 2.197849
```

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/nhejazi/txshift/issues).

-----

## Funding

The development of this software was supported in part through a grant
from the National Library of Medicine of the National Institutes of
Health (T32 LM012417).

-----

## License

© 2017-2018 [Nima S. Hejazi](http://nimahejazi.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    
    Copyright (c) 2017-2018 Nima S. Hejazi
    
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

-----

## References

<div id="refs" class="references">

<div id="ref-diaz2017stochastic">

Díaz, Iván, and Mark J van der Laan. 2017. “Stochastic Treatment
Regimes.” In *Targeted Learning in Data Science: Causal Inference for
Complex Longitudinal Studies*, 167–80. Springer Science & Business
Media.

</div>

<div id="ref-munoz2012population">

Muñoz, Iván Díaz, and Mark J van der Laan. 2012. “Population
Intervention Causal Effects Based on Stochastic Interventions.”
*Biometrics* 68 (2). Wiley Online Library:541–49.

</div>

<div id="ref-vdl2011targeted">

van der Laan, Mark J, and Sherri Rose. 2011. *Targeted Learning: Causal
Inference for Observational and Experimental Data*. Springer Science &
Business Media.

</div>

<div id="ref-vdl2017targeted">

———. 2017. *Targeted Learning in Data Science: Causal Inference for
Complex Longitudinal Studies*. Springer Science & Business Media.

</div>

</div>
