
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`txshift`

[![Travis-CI Build
Status](https://travis-ci.org/nhejazi/txshift.svg?branch=master)](https://travis-ci.org/nhejazi/txshift)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/nhejazi/txshift?branch=master&svg=true)](https://ci.appveyor.com/project/nhejazi/txshift)
[![Coverage
Status](https://img.shields.io/codecov/c/github/nhejazi/txshift/master.svg)](https://codecov.io/github/nhejazi/txshift?branch=master)
[![CRAN](http://www.r-pkg.org/badges/version/txshift)](http://www.r-pkg.org/pkg/txshift)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/txshift)](https://CRAN.R-project.org/package=txshift)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> Targeted Learning of the Causal Effects of Stochastic Interventions

**Authors:** [Nima Hejazi](https://nimahejazi.org) and [David
Benkeser](https://www.benkeserstatistics.com/)

-----

## What’s `txshift`?

The `txshift` R package is designed to compute targeted maximum
likelihood (TML) estimates of the counterfactual mean of an outcome
under stochastic mechanisms for treatment assignment and related causal
parameters (Díaz and van der Laan (2012)). In particular, `txshift`
implements and builds upon a simplified algorithm for the TML estimator
proposed by Díaz and van der Laan (2018).

For many practical applications (e.g., vaccine efficacy trials), it is
often the case that the observed data structure is generated under a
two-phase sampling mechanism (i.e., a two-stage design). In such cases,
TML estimators must be augmented to exhibit efficiency in spite of the
challenges induced by the censoring process. An appropriate augmentation
procedure was first proposed by Rose and van der Laan (2011), who
proposed the use of inverse probability of censoring weights (IPCW)
alongside an augmentation of the relevant efficient influence function.
`txshift` extends this approach to computing one-step and IPCW-TML
estimators of the counterfactual mean under a stochastic treatment
regime.

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
library(txshift)
library(condensier)
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
Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))

# fit the TMLE
tmle_shift <- tmle_txshift(W = W, A = A, Y = Y, delta = 0.5,
                           g_fit_args = list(fit_type = "glm",
                                             nbins = 25,
                                             bin_method = "dhist",
                                             bin_estimator =
                                               speedglmR6$new(),
                                             parfit = FALSE),
                           Q_fit_args = list(fit_type = "glm",
                                             glm_formula = "Y ~ .")
                          )

# conveniently summarize the results
summary(tmle_shift)
#>        lwr_ci     param_est        upr_ci     param_var      eif_mean 
#>      0.722106      0.752212       0.78005      0.000219 -3.381496e-17 
#>        n_iter 
#>             0

# now, let's introduce a censoring process (for two-stage sampling)
C <- rbinom(n_obs, 1, plogis(W + Y))

# fit an IPCW-TMLE to account for this censoring process:
ipcwtmle_shift <- tmle_txshift(W = W, A = A, Y = Y, delta = 0.5,
                               C = C, V = c("W", "Y"),
                               max_iter = 10,  # limit iterations for speed
                               ipcw_fit_args = list(fit_type = "glm"),
                               g_fit_args = list(fit_type = "glm",
                                                 nbins = 25,
                                                 bin_method = "dhist",
                                                 bin_estimator =
                                                   speedglmR6$new(),
                                                 parfit = FALSE),
                               Q_fit_args = list(fit_type = "glm",
                                                 glm_formula = "Y ~ ."),
                               eif_reg_type = "glm"
                              )

# conveniently summarize the results
summary(ipcwtmle_shift)
#>        lwr_ci     param_est        upr_ci     param_var      eif_mean 
#>      0.728859      0.764822      0.797341      0.000306 -2.964196e-04 
#>        n_iter 
#>             1
```

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/nhejazi/txshift/issues).

-----

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/nhejazi/txshift/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

-----

## Citation

After using the `txshift` R package, please cite the following:

``` 
    @manual{hejazi2019txshift,
      author = {Hejazi, Nima S and Benkeser, David C},
      title = {txshift: {Targeted Learning} of the Causal Effects of
        Stochastic Interventions in {R}},
      year  = {2019},
      url = {https://github.com/nhejazi/txshift},
      note = {R package version 0.2.2}
    }
```

-----

## Related

  - [R/`tmle3shift`](https://github.com/tlverse/tmle3shift) - An R
    package providing an independent implementation of the same core
    routines for the TML estimation procedure and statistical
    methodology as is made available here, through reliance on a unified
    interface for Targeted Learning provided by the
    [`tmle3`](https://github.com/tlverse/tmle3) engine of the [`tlverse`
    ecosystem](https://github.com/tlverse).

  - [R/`medshift`](https://github.com/nhejazi/medshift) - An R package
    providing facilities to estimate the causal effect of stochastic
    treatment regimes in the mediation setting, including classical
    (IPW) and augmented double robust (one-step) estimators. This is an
    implementation of the methodology explored in Díaz and Hejazi
    (2019).

  - [R/`condensier`](https://github.com/osofr/condensier) - Estimation
    of the treatment mechanism component of this parameter requires
    conditional density estimation, which is implemented rather
    generally (for a variety of hazard regression strategies and
    regression functions) in this package. This package implements a
    variant of the methodology proposed in Díaz and van der Laan (2011).

  - [R/`haldensify`](https://github.com/nhejazi/haldensify) - A minimal
    package for estimating the conditional density treatment mechanism
    component of this parameter based on using the [highly adaptive
    lasso](https://github.com/tlverse/hal9001) for the required pooled
    hazard regression. This package implements the methdology proposed
    in Díaz and van der Laan (2011).

-----

## Funding

The development of this software was supported in part through a grant
from the National Institutes of Health: [T32
LM012417-02](https://projectreporter.nih.gov/project_info_description.cfm?aid=9248418&icde=37849831&ddparam=&ddvalue=&ddsub=&cr=1&csb=default&cs=ASC&pball=).

-----

## License

© 2017-2019 [Nima S. Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    
    Copyright (c) 2017-2019 Nima S. Hejazi
    
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

<div id="ref-diaz2019causal">

Díaz, Iván, and Nima S Hejazi. 2019. “Causal Mediation Analysis for
Stochastic Interventions.” *Submitted*.
<https://arxiv.org/abs/1901.02776>.

</div>

<div id="ref-diaz2011super">

Díaz, Iván, and Mark J van der Laan. 2011. “Super Learner Based
Conditional Density Estimation with Application to Marginal Structural
Models.” *The International Journal of Biostatistics* 7 (1). De Gruyter:
1–20.

</div>

<div id="ref-diaz2012population">

———. 2012. “Population Intervention Causal Effects Based on Stochastic
Interventions.” *Biometrics* 68 (2). Wiley Online Library: 541–49.

</div>

<div id="ref-diaz2018stochastic">

———. 2018. “Stochastic Treatment Regimes.” In *Targeted Learning in Data
Science: Causal Inference for Complex Longitudinal Studies*, 167–80.
Springer Science & Business Media.

</div>

<div id="ref-rose2011targeted2sd">

Rose, Sherri, and Mark J van der Laan. 2011. “A Targeted Maximum
Likelihood Estimator for Two-Stage Designs.” *The International Journal
of Biostatistics* 7 (1): 1–21.

</div>

</div>
