
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`txshift`

<!-- badges: start -->

[![R-CMD-check](https://github.com/nhejazi/txshift/workflows/R-CMD-check/badge.svg)](https://github.com/nhejazi/txshift/actions)
[![Coverage
Status](https://img.shields.io/codecov/c/github/nhejazi/txshift/master.svg)](https://codecov.io/github/nhejazi/txshift?branch=master)
[![CRAN](https://www.r-pkg.org/badges/version/txshift)](https://www.r-pkg.org/pkg/txshift)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/txshift)](https://CRAN.R-project.org/package=txshift)
[![CRAN total
downloads](http://cranlogs.r-pkg.org/badges/grand-total/txshift)](https://CRAN.R-project.org/package=txshift)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT
license](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4070042.svg)](https://doi.org/10.5281/zenodo.4070042)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02447/status.svg)](https://doi.org/10.21105/joss.02447)
<!-- badges: end -->

> Efficient Estimation of the Causal Effects of Stochastic Interventions

**Authors:** [Nima Hejazi](https://nimahejazi.org) and [David
Benkeser](https://www.sph.emory.edu/faculty/profile/#!dbenkes)

-----

## What’s `txshift`?

The `txshift` R package is designed to provide facilities for the
construction of efficient estimators of a causal parameter defined as
the counterfactual mean of an outcome under stochastic mechanisms for
treatment assignment (Dı́az and van der Laan 2012). `txshift` implements
and builds upon a simplified algorithm for the targeted maximum
likelihood (TML) estimator of such a causal parameter, originally
proposed by Dı́az and van der Laan (2018), and makes use of analogous
machinery to compute an efficient one-step estimator (Pfanzagl and
Wefelmeyer 1985). `txshift` integrates with the [`sl3`
package](https://github.com/tlverse/sl3) (Coyle et al. 2020) to allow
for ensemble machine learning to be leveraged in the estimation
procedure.

For many practical applications (e.g., vaccine efficacy trials),
observed data is often subject to a two-phase sampling mechanism (i.e.,
through the use of a two-stage design). In such cases, efficient
estimators (of both varieties) must be augmented to construct unbiased
estimates of the population-level causal parameter. Rose and van der
Laan (2011) first introduced an augmentation procedure that relies on
introducing inverse probability of censoring (IPC) weights directly to
an appropriate loss function or to the efficient influence function
estimating equation. `txshift` extends this approach to compute
IPC-weighted one-step and TML estimators of the counterfactual mean
outcome under a shift stochastic treatment regime. The package is
designed to implement the statistical methodology described in Hejazi et
al. (2020) and extensions thereof.

-----

## Installation

For standard use, we recommend installing the package from
[CRAN](https://CRAN.R-project.org/package=txshift) via

``` r
install.packages("txshift")
```

*Note:* If `txshift` is installed from
[CRAN](https://CRAN.R-project.org/package=txshift), the `sl3`, an
enhancing dependency that allows ensemble machine learning to be used
for nuisance parameter estimation, won’t be included. We highly
recommend additionally installing `sl3` from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("tlverse/sl3@master")
```

For the latest features, install the most recent *stable version* of
`txshift` from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("nhejazi/txshift@master")
```

To contribute, install the *development version* of `txshift` from
GitHub via [`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("nhejazi/txshift@devel")
```

-----

## Example

To illustrate how `txshift` may be used to ascertain the effect of a
treatment, consider the following example:

``` r
library(txshift)
#> txshift v0.3.6: Efficient Estimation of the Causal Effects of Stochastic
#> Interventions
library(sl3)
set.seed(429153)

# simulate simple data
n_obs <- 500
W <- replicate(2, rbinom(n_obs, 1, 0.5))
A <- rnorm(n_obs, mean = 2 * W, sd = 1)
Y <- rbinom(n_obs, 1, plogis(A + W + rnorm(n_obs, mean = 0, sd = 1)))

# now, let's introduce a a two-stage sampling process
C_samp <- rbinom(n_obs, 1, plogis(W + Y))

# fit the full-data TMLE (ignoring two-phase sampling)
tmle <- txshift(
  W = W, A = A, Y = Y, delta = 0.5,
  estimator = "tmle",
  g_exp_fit_args = list(
    fit_type = "sl",
    sl_learners_density = Lrnr_density_hse$new(Lrnr_hal9001$new())
  ),
  Q_fit_args = list(fit_type = "glm", glm_formula = "Y ~ .")
)
tmle
#> Counterfactual Mean of Shifted Treatment
#> Intervention: Treatment + 0.5
#> txshift Estimator: tmle
#> Estimate: 0.7685
#> Std. Error: 0.019
#> 95% CI: [0.7292, 0.8037]
# fit a full-data one-step estimator for comparison (again, no sampling)
os <- txshift(
  W = W, A = A, Y = Y, delta = 0.5,
  estimator = "onestep",
  g_exp_fit_args = list(
    fit_type = "sl",
    sl_learners_density = Lrnr_density_hse$new(Lrnr_hal9001$new())
  ),
  Q_fit_args = list(fit_type = "glm", glm_formula = "Y ~ .")
)
os
#> Counterfactual Mean of Shifted Treatment
#> Intervention: Treatment + 0.5
#> txshift Estimator: onestep
#> Estimate: 0.7685
#> Std. Error: 0.019
#> 95% CI: [0.7292, 0.8037]
# fit an IPCW-TMLE to account for the two-phase sampling process
tmle_ipcw <- txshift(
  W = W, A = A, Y = Y, delta = 0.5, C_samp = C_samp, V = c("W", "Y"),
  estimator = "tmle", max_iter = 5, eif_reg_type = "glm",
  samp_fit_args = list(fit_type = "glm"),
  g_exp_fit_args = list(
    fit_type = "sl",
    sl_learners_density = Lrnr_density_hse$new(Lrnr_hal9001$new())
  ),
  Q_fit_args = list(fit_type = "glm", glm_formula = "Y ~ .")
)
tmle_ipcw
#> Counterfactual Mean of Shifted Treatment
#> Intervention: Treatment + 0.5
#> txshift Estimator: tmle
#> Estimate: 0.7603
#> Std. Error: 0.0204
#> 95% CI: [0.718, 0.798]
# compare with an IPCW-agumented one-step estimator under two-phase sampling
os_ipcw <- txshift(
  W = W, A = A, Y = Y, delta = 0.5, C_samp = C_samp, V = c("W", "Y"),
  estimator = "onestep", eif_reg_type = "glm",
  samp_fit_args = list(fit_type = "glm"),
  g_exp_fit_args = list(
    fit_type = "sl",
    sl_learners_density = Lrnr_density_hse$new(Lrnr_hal9001$new())
  ),
  Q_fit_args = list(fit_type = "glm", glm_formula = "Y ~ .")
)
os_ipcw
#> Counterfactual Mean of Shifted Treatment
#> Intervention: Treatment + 0.5
#> txshift Estimator: onestep
#> Estimate: 0.7601
#> Std. Error: 0.0204
#> 95% CI: [0.7178, 0.7979]
```

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/nhejazi/txshift/issues). Further
details on filing issues are provided in our [contribution
guidelines](https://github.com/nhejazi/txshift/blob/master/CONTRIBUTING.md).

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
    @article{hejazi2020efficient,
      author = {Hejazi, Nima S and {van der Laan}, Mark J and Janes, Holly
        E and Gilbert, Peter B and Benkeser, David C},
      title = {Efficient nonparametric inference on the effects of
        stochastic interventions under two-phase sampling, with
        applications to vaccine efficacy trials},
      year = {2020},
      doi = {10.1111/biom.13375},
      url = {https://doi.org/10.1111/biom.13375},
      journal = {Biometrics},
      publisher = {Wiley Online Library}
    }

    @article{hejazi2020txshift-joss,
      author = {Hejazi, Nima S and Benkeser, David C},
      title = {{txshift}: Efficient estimation of the causal effects of
        stochastic interventions in {R}},
      year  = {2020},
      doi = {10.21105/joss.02447},
      url = {https://10.21105.joss.02447},
      journal = {Journal of Open Source Software},
      publisher = {The Open Journal}
    }

    @software{hejazi2020txshift-rpkg,
      author = {Hejazi, Nima S and Benkeser, David C},
      title = {{txshift}: Efficient Estimation of the Causal Effects of
        Stochastic Interventions},
      year  = {2020},
      doi = {10.5281/zenodo.4070042},
      url = {https://CRAN.R-project.org/package=txshift},
      note = {R package version 0.3.4}
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
    implementation of the methodology explored by Dı́az and Hejazi
    (2020).

  - [R/`haldensify`](https://github.com/nhejazi/haldensify) - A minimal
    package for estimating the conditional density treatment mechanism
    component of this parameter based on using the [highly adaptive
    lasso](https://github.com/tlverse/hal9001) (Coyle, Hejazi, and van
    der Laan 2020; Hejazi, Coyle, and van der Laan 2020) in combination
    with a pooled hazard regression. This package implements a variant
    of the approach advocated by Dı́az and van der Laan (2011).

-----

## Funding

The development of this software was supported in part through grants
from the National Library of Medicine (award number [T32
LM012417](https://reporter.nih.gov/project-details/9248418)) and the
National Institute of Allergy and Infectious Diseases (award number [R01
AI074345](https://reporter.nih.gov/project-details/9926564)) of the
National Institutess of Health.

-----

## License

© 2017-2021 [Nima S. Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    
    Copyright (c) 2017-2021 Nima S. Hejazi
    
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

<div id="ref-coyle2020sl3">

Coyle, Jeremy R, Nima S Hejazi, Ivana Malenica, and Oleg Sofrygin. 2020.
*sl3: Modern Pipelines for Machine Learning and Super Learning*.
<https://github.com/tlverse/sl3>.
<https://doi.org/10.5281/zenodo.1342293>.

</div>

<div id="ref-coyle2020hal9001-rpkg">

Coyle, Jeremy R, Nima S Hejazi, and Mark J van der Laan. 2020. *hal9001:
The Scalable Highly Adaptive Lasso*.
<https://github.com/tlverse/hal9001>.
<https://doi.org/10.5281/zenodo.3558313>.

</div>

<div id="ref-diaz2020causal">

Dı́az, Iván, and Nima S Hejazi. 2020. “Causal Mediation Analysis for
Stochastic Interventions.” *Journal of the Royal Statistical Society:
Series B (Statistical Methodology)* 82 (3): 661–83.
<https://doi.org/10.1111/rssb.12362>.

</div>

<div id="ref-diaz2011super">

Dı́az, Iván, and Mark J van der Laan. 2011. “Super Learner Based
Conditional Density Estimation with Application to Marginal Structural
Models.” *The International Journal of Biostatistics* 7 (1): 1–20.

</div>

<div id="ref-diaz2012population">

———. 2012. “Population Intervention Causal Effects Based on Stochastic
Interventions.” *Biometrics* 68 (2): 541–49.

</div>

<div id="ref-diaz2018stochastic">

———. 2018. “Stochastic Treatment Regimes.” In *Targeted Learning in Data
Science: Causal Inference for Complex Longitudinal Studies*, 167–80.
Springer Science & Business Media.

</div>

<div id="ref-hejazi2020hal9001-joss">

Hejazi, Nima S, Jeremy R Coyle, and Mark J van der Laan. 2020. “hal9001:
Scalable Highly Adaptive Lasso Regression in R.” *Journal of Open Source
Software* 5 (53): 2526. <https://doi.org/10.21105/joss.02526>.

</div>

<div id="ref-hejazi2020efficient">

Hejazi, Nima S, Mark J van der Laan, Holly E Janes, Peter B Gilbert, and
David C Benkeser. 2020. “Efficient Nonparametric Inference on the
Effects of Stochastic Interventions Under Two-Phase Sampling, with
Applications to Vaccine Efficacy Trials.” *Biometrics*.
<https://doi.org/10.1111/biom.13375>.

</div>

<div id="ref-pfanzagl1985contributions">

Pfanzagl, J, and W Wefelmeyer. 1985. “Contributions to a General
Asymptotic Statistical Theory.” *Statistics & Risk Modeling* 3 (3-4):
379–88.

</div>

<div id="ref-rose2011targeted2sd">

Rose, Sherri, and Mark J van der Laan. 2011. “A Targeted Maximum
Likelihood Estimator for Two-Stage Designs.” *The International Journal
of Biostatistics* 7 (1): 1–21.

</div>

</div>
