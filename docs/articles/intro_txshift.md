# Evaluating Causal Effects of Modified Treatment Policies

## Introduction

Stochastic treatment regimes constitute a flexible framework for
evaluating the effects of continuous-valued exposures/treatments.
*Modified treatment policies*, one such technique within this framework,
examine the effects attributable to shifting the observed (“natural”)
value of a treatment, usually up or down by some scalar $`\delta`$. The
`txshift` package implements algorithms for computing one-step or
targeted minimum loss-based (TML) estimates of the counterfactual means
induced by additive modified treatment policies (MTPs), defined by a
shifting function $`\delta(A,W)`$. For a technical presentation, the
interested reader is invited to consult Dı́az and van der Laan (2018) or
the earlier work of Dı́az and van der Laan (2012) and Haneuse and
Rotnitzky (2013). For background on Targeted Learning, consider
consulting van der Laan and Rose (2011), van der Laan and Rose (2018),
and van der Laan et al. (2022).

To start, let’s load the packages we’ll need and set a seed for
simulation:

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`data.table`](https://r-datatable.com)`)`

    ## 
    ## Attaching package: 'data.table'

    ## The following object is masked from 'package:base':
    ## 
    ##     %notin%

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`haldensify`](https://codex.nimahejazi.org/haldensify/)`)`

    ## haldensify v0.2.8: Highly Adaptive Lasso Conditional Density Estimation

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`txshift`](https://github.com/nhejazi/txshift)`)`

    ## txshift v0.4.0: Efficient Estimation of the Causal Effects of Stochastic Interventions

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``11249``)`

------------------------------------------------------------------------

## Data and Notation

We’ll consider $`n`$ observed units $`O_1, \ldots, O_n`$, where each
random variable $`O = (W, A, Y)`$ corresponds to the data available on a
single unit. Within $`O`$, $`W`$ denotes baseline covariates (e.g., age,
biological sex, BMI), $`A \in
\mathbb{R}`$ a continuous-valued exposure (e.g., dosage of nutritional
supplements taken), and $`Y`$ an outcome of interest (e.g., disease
status). To minimize unjustifiable assumptions, we let
$`O \sim \mathcal{P} \in \mathcal{M}`$, where $`\mathcal{P}`$ is simply
any distribution within the nonparametric statistical model
$`\mathcal{M}`$. To formalize the definition of stochastic interventions
and their corresponding causal effects, we consider a nonparametric
structural equation model (NPSEM), introduced by Pearl (2000), to define
how the system changes under interventions of interest:
``` math
\begin{align*}\label{eqn:npsem}
  W &= f_W(U_W) \\ A &= f_A(W, U_A) \\ Y &= f_Y(A, W, U_Y),
\end{align*}
```
We denote the observed data structure $`O = (W, A, Y)`$

Assuming that the distribution of $`A`$ conditional on $`W = w`$ has
support in the interval $`(l(w), u(w))`$ – for convenience, we assume
that the minimum natural value of treatment $`A`$ for an individual with
covariates $`W = w`$ is $`l(w)`$, while, similarly, the maximum is
$`u(w)`$ – a simple MTP based on a shift $`\delta`$, is
``` math
\begin{equation}\label{eqn:shift}
  \delta(a, w) =
  \begin{cases}
    a - \delta & \text{if } a > l(w) + \delta \\
    a & \text{if } a \leq l(w) + \delta,
  \end{cases}
\end{equation}
```
where $`0 \leq \delta`$ is an arbitrary pre-specified value that defines
the degree to which the observed value $`A`$ is to be shifted, where
possible.

In case-cohort studies, it is common practice to make use of
outcome-dependent two-phase sampling designs, which allow for expensive
measurements made on the exposure (e.g., genomic sequencing of immune
markers) to be avoided. As a complication, such sampling schemes alter
the observed data structure from the simpler $`O = (W, A, Y)`$ to
$`O = (W, \Delta A, Y, \Delta)`$, where the sampling indicator
$`\Delta`$ may itself be a function of the variables $`\{W, Y\}`$. In
this revised data structure, the value of $`A`$ is only observed for
units in the two-phase sample, for whom $`\Delta = 1`$. Hejazi et al.
(2020) provide a detailed investigation of the methodological details of
efficient estimation under such designs in the context of vaccine
efficacy trials; their work was used in the analysis of immune
correlates of protection for HIV-1 (Hejazi et al. 2020) and COVID-19
(Gilbert et al. 2021). Of course, one may also account for loss to
follow-up (i.e., censoring), which the `txshift` package supports
(through the `C_cens` argument of the eponymous `txshift` function),
though we avoid this complication in our subsequent examples in the
interest of clarity of exposition. Corrections for both censoring and
two-phase sampling make use of inverse probability of censoring
weighting (IPCW), leading to IPCW-augmented one-step and TML estimators.

### Simulate Data

\
`# parameters for simulation example`\
`n_obs`` ``<-`` ``200`` ``# number of observations`\
\
`# baseline covariate -- simple, binary`\
`W`` ``<-`` `[`rbinom`](https://rdrr.io/r/stats/Binomial.html)`(``n_obs``, ``1``, ``0.5``)`\
\
`# create treatment based on baseline W`\
`A`` ``<-`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_obs``, mean ``=`` ``2`` ``*`` ``W``, sd ``=`` ``0.5``)`\
\
`# create outcome as a linear function of A, W + white noise`\
`Y`` ``<-`` ``A`` ``+`` ``W`` ``+`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n_obs``, mean ``=`` ``0``, sd ``=`` ``1``)`\
\
`# two-phase sampling based on covariates`\
`Delta_samp`` ``<-`` `[`rbinom`](https://rdrr.io/r/stats/Binomial.html)`(``n_obs``, ``1``, `[`plogis`](https://rdrr.io/r/stats/Logistic.html)`(``W``)``)`\
\
`# treatment shift parameter`\
`delta`` ``<-`` ``0.5`

## Estimating the Effects of Additive MTPs

The simplest way to compute an efficient estimator for an additive MTP
is to fit each of the nuisance parameters internally. This procedure can
be sped up by using generalized linear models (GLMs) to fit the outcome
regression $`Q_n`$. The
[`txshift()`](https://codex.nimahejazi.org/txshift/reference/txshift.md)
function provides a simple interface.

\
`est_shift`` ``<-`` `[`txshift`](https://codex.nimahejazi.org/txshift/reference/txshift.md)`(`\
`  W ``=`` ``W``, A ``=`` ``A``, Y ``=`` ``Y``, delta ``=`` ``delta``,`\
`  g_exp_fit_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`    fit_type ``=`` ``"hal"``, n_bins ``=`` ``5``,`\
`    grid_type ``=`` ``"equal_mass"``,`\
`    lambda_seq ``=`` `[`exp`](https://rdrr.io/r/base/Log.html)`(`[`seq`](https://rdrr.io/r/base/seq.html)`(``-``1``, ``-``10``, length ``=`` ``100``)``)`\
`  ``)``,`\
`  Q_fit_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`    fit_type ``=`` ``"glm"``,`\
`    glm_formula ``=`` ``"Y ~ .^2"`\
`  ``)`\
`)`\
`est_shift`

    ## Counterfactual Mean of Shifted Treatment

    ## Intervention: Treatment + 0.5

    ## txshift Estimator: tmle

    ## Estimate: 1.8623

    ## Std. Error: 0.1542

    ## 95% CI: [1.56, 2.1646]

### Interlude: Ensemble Machine Learning with `sl3`

To easily incorporate ensemble machine learning into the estimation
procedure, `txshift` integrates with the [`sl3` R
package](https://tlverse.org/sl3/) (Coyle, Hejazi, Malenica, et al.
2022). For a complete guide on using the `sl3` R package, consider
consulting [the chapter on Super
Learning](https://tlverse.org/tlverse-handbook/sl3.html) in van der Laan
et al. (2022).

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`sl3`](https://tlverse.org/sl3)`)`\
\
`# SL learners to be used for most fits (e.g., IPCW, outcome regression)`\
`mean_learner`` ``<-`` `[`Lrnr_mean`](https://tlverse.org/sl3/reference/Lrnr_mean.html)`$``new``(``)`\
`glm_learner`` ``<-`` `[`Lrnr_glm`](https://tlverse.org/sl3/reference/Lrnr_glm.html)`$``new``(``)`\
`rf_learner`` ``<-`` `[`Lrnr_ranger`](https://tlverse.org/sl3/reference/Lrnr_ranger.html)`$``new``(``)`\
`Q_lib`` ``<-`` `[`Stack`](https://tlverse.org/sl3/reference/Stack.html)`$``new``(``mean_learner``, ``glm_learner``, ``rf_learner``)`\
`sl_learner`` ``<-`` `[`Lrnr_sl`](https://tlverse.org/sl3/reference/Lrnr_sl.html)`$``new``(``learners ``=`` ``Q_lib``, metalearner ``=`` `[`Lrnr_nnls`](https://tlverse.org/sl3/reference/Lrnr_nnls.html)`$``new``(``)``)`\
\
`# SL learners for fitting the generalized propensity score fit`\
`hose_learner`` ``<-`` `[`make_learner`](https://tlverse.org/sl3/reference/Lrnr_base.html)`(``Lrnr_density_semiparametric``,`\
`  mean_learner ``=`` ``glm_learner`\
`)`\
`hese_learner`` ``<-`` `[`make_learner`](https://tlverse.org/sl3/reference/Lrnr_base.html)`(``Lrnr_density_semiparametric``,`\
`  mean_learner ``=`` ``rf_learner``,`\
`  var_learner ``=`` ``glm_learner`\
`)`\
`g_lib`` ``<-`` `[`Stack`](https://tlverse.org/sl3/reference/Stack.html)`$``new``(``hose_learner``, ``hese_learner``)`\
`sl_learner_density`` ``<-`` `[`Lrnr_sl`](https://tlverse.org/sl3/reference/Lrnr_sl.html)`$``new``(`\
`  learners ``=`` ``g_lib``,`\
`  metalearner ``=`` `[`Lrnr_solnp_density`](https://tlverse.org/sl3/reference/Lrnr_solnp_density.html)`$``new``(``)`\
`)`

### Efficient Effect Estimates with Machine Learning

Using the framework provided by the [`sl3`
package](https://github.com/tlverse/sl3/), the nuisance functions
required for our efficient estimators may be fit with ensemble machine
learning. The Super Learner algorithm (van der Laan et al. 2007)
implemented in `sl3` uses the asymptotic optimality of V-fold
cross-validation (Dudoit and van der Laan 2005; van der Laan et al.
2004; van der Vaart et al. 2006) to select an optimal prediction
functions from a library or to construct an optimal combination of
prediction functions.

\
`est_shift_sl`` ``<-`` `[`txshift`](https://codex.nimahejazi.org/txshift/reference/txshift.md)`(`\
`  W ``=`` ``W``, A ``=`` ``A``, Y ``=`` ``Y``, delta ``=`` ``delta``,`\
`  g_exp_fit_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`    fit_type ``=`` ``"sl"``,`\
`    sl_learners_density ``=`` ``sl_learner_density`\
`  ``)``,`\
`  Q_fit_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`    fit_type ``=`` ``"sl"``,`\
`    sl_learners ``=`` ``sl_learner`\
`  ``)`\
`)`\
`est_shift_sl`

## Estimating the Effects of Additive MTPs Under Two-Phase Sampling

In case-cohort studies in which two-phase sampling is used, the data
structure takes the from $`O = (W, \Delta A, Y, \Delta)`$ as previously
discussed. Under such sampling, the
[`txshift()`](https://codex.nimahejazi.org/txshift/reference/txshift.md)
function may still be used to estimate the causal effect of an additive
MTP – only a few additional arguments need to be specified:

\
`est_shift_ipcw`` ``<-`` `[`txshift`](https://codex.nimahejazi.org/txshift/reference/txshift.md)`(`\
`  W ``=`` ``W``, A ``=`` ``A``, Y ``=`` ``Y``, delta ``=`` ``delta``,`\
`  C_samp ``=`` ``Delta_samp``, V ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"W"``, ``"Y"``)``,`\
`  samp_fit_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``fit_type ``=`` ``"glm"``)``,`\
`  g_exp_fit_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`    fit_type ``=`` ``"hal"``, n_bins ``=`` ``5``,`\
`    grid_type ``=`` ``"equal_mass"``,`\
`    lambda_seq ``=`` `[`exp`](https://rdrr.io/r/base/Log.html)`(`[`seq`](https://rdrr.io/r/base/seq.html)`(``-``1``, ``-``10``, length ``=`` ``100``)``)`\
`  ``)``,`\
`  Q_fit_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`    fit_type ``=`` ``"glm"``,`\
`    glm_formula ``=`` ``"Y ~ .^2"`\
`  ``)``,`\
`  eif_reg_type ``=`` ``"glm"`\
`)`\
`est_shift_ipcw`

    ## Counterfactual Mean of Shifted Treatment

    ## Intervention: Treatment + 0.5

    ## txshift Estimator: tmle

    ## Estimate: 1.958

    ## Std. Error: 0.1396

    ## 95% CI: [1.6844, 2.2316]

Note that we specify a few additional arguments in the call to
[`txshift()`](https://codex.nimahejazi.org/txshift/reference/txshift.md),
including `C_samp`, the indicator of inclusion in the two-phase sample;
`V`, the set of other variables that may affect the sampling decision
(in this case, both the baseline covariates and the outcome);
`samp_fit_args`, which indicates how the sampling mechanism ought to be
estimated; and `eif_reg_type`, which indicates how a particular
reduced-dimension nuisance regression ought to be estimated (see Rose
and van der Laan (2011) and Hejazi et al. (2020) for details). This last
argument only has options for using a GLM or the highly adaptive lasso
(HAL), a nonparametric regression estimator, using the [`hal9001`
package](https://github.com/tlverse/hal9001/) \[Coyle, Hejazi, Phillips,
et al. (2022); hejazi2020hal9001-joss\]. In practice, we recommend
leaving this argument to the default and fitting this nuisance function
with HAL; however, this is significantly more costly computationally.

## Statistical Inference for One-step and TML Estimators

The efficient estimators implemented in `txshift` are asymptotically
linear; thus, the estimator $`\psi_n`$ converges to the true parameter
value $`\psi_0`$:
``` math
\psi_n - \psi_0 = (P_n - P_0) \cdot D(\bar{Q}_n^{\star}, g_n) +
  R(\hat{P}^{\star}, P_0),
```
which yields
``` math
\psi_n - \psi_0 = (P_n - P_0) \cdot D(P_0) + o_P \left( \frac{1}{\sqrt{n}}
 \right),
```
provided the following conditions,

1.  if $`D(\bar{Q}_n^{\star}, g_n)`$ converges to $`D(P_0)`$ in
    $`L_2(P_0)`$ norm;
2.  the size of the class of functions considered for estimation of
    $`\bar{Q}_n^{\star}`$ and $`g_n`$ is bounded (technically,
    $`\exists \mathcal{F}`$ such that
    $`D(\bar{Q}_n^{\star}, g_n) \in \mathcal{F}`$ with high probability,
    where $`\mathcal{F}`$ is a Donsker class); and
3.  the remainder term $`R(\hat{P}^{\star}, P_0)`$ decays as $`o_P
    \left( \frac{1}{\sqrt{n}} \right)`$.

By the central limit theorem, the estimators then have a Gaussian
limiting distribution,
``` math
\sqrt{n}(\psi_n - \psi) \to N(0, V(D(P_0))),
```
where $`V(D(P_0))`$ is the variance of the efficient influence function
(or canonical gradient).

The above implies that $`\psi_n`$ is a $`\sqrt{n}`$-consistent estimator
of $`\psi`$, that it is asymptotically normal (as given above), and that
it is locally efficient. This allows us to build Wald-type confidence
intervals in a straightforward manner:

``` math
\psi_n \pm z_{\alpha} \cdot \frac{\sigma_n}{\sqrt{n}},
```
where $`\sigma_n^2`$ is an estimator of $`V(D(P_0))`$. The estimator
$`\sigma_n^2`$ may be obtained using the bootstrap or computed directly
via the following

``` math
\sigma_n^2 = \frac{1}{n} \sum_{i = 1}^{n} D^2(\bar{Q}_n^{\star}, g_n)(O_i).
```

Such confidence intervals may easily be created with the `confint`
method:

\
`(``ci_est_shift`` ``<-`` `[`confint`](https://rdrr.io/r/stats/confint.html)`(``est_shift``)``)`

    ##   lwr_ci      est   upr_ci 
    ## 1.559997 1.862286 2.164576

## *Advanced Usage:* User-Specified Nuisance Regressions

In some special cases it may be useful for the experienced user to
compute the treatment mechanism, censoring mechanism, outcome
regression, and sampling mechanism fits separately (i.e., outside of the
`txshift` wrapper function), instead applying the wrapper only to
construct an efficient one-step or TML estimator. In such cases, the
optional arguments ending in `_ext`.

\
`# compute censoring mechanism and produce IPC weights externally`\
`pi_mech`` ``<-`` `[`plogis`](https://rdrr.io/r/stats/Logistic.html)`(``W``)`\
`ipcw_out`` ``<-`` ``pi_mech`\
\
`# compute treatment mechanism (propensity score) externally`\
`## first, produce the down-shifted treatment data`\
`gn_downshift`` ``<-`` `[`dnorm`](https://rdrr.io/r/stats/Normal.html)`(``A`` ``-`` ``delta``, mean ``=`` ``tx_mult`` ``*`` ``W``, sd ``=`` ``1``)`\
`## next, initialize and produce the up-shifted treatment data`\
`gn_upshift`` ``<-`` `[`dnorm`](https://rdrr.io/r/stats/Normal.html)`(``A`` ``+`` ``delta``, mean ``=`` ``tx_mult`` ``*`` ``W``, sd ``=`` ``1``)`\
`## now, initialize and produce the up-up-shifted (2 * delta) treatment data`\
`gn_upupshift`` ``<-`` `[`dnorm`](https://rdrr.io/r/stats/Normal.html)`(``A`` ``+`` ``2`` ``*`` ``delta``, mean ``=`` ``tx_mult`` ``*`` ``W``, sd ``=`` ``1``)`\
`## then, initialize and produce the un-shifted treatment data`\
`gn_noshift`` ``<-`` `[`dnorm`](https://rdrr.io/r/stats/Normal.html)`(``A``, mean ``=`` ``tx_mult`` ``*`` ``W``, sd ``=`` ``1``)`\
`## finally, put it all together into an object like what's produced internally`\
`gn_out`` ``<-`` `[`as.data.table`](https://rdrr.io/pkg/data.table/man/as.data.table.html)`(`[`cbind`](https://rdrr.io/r/base/cbind.html)`(`\
`  ``gn_downshift``, ``gn_noshift``, ``gn_upshift``,`\
`  ``gn_upupshift`\
`)``)``[``C`` ``==`` ``1``, ``]`\
[`setnames`](https://rdrr.io/pkg/data.table/man/setattr.html)`(``gn_out``, `[`c`](https://rdrr.io/r/base/c.html)`(``"downshift"``, ``"noshift"``, ``"upshift"``, ``"upupshift"``)``)`\
\
`# compute outcome regression externally`\
`# NOTE: transform Y to lie in the unit interval and bound predictions such that`\
`#       no values fall near the bounds of the interval`\
`Qn_noshift`` ``<-`` ``(``W`` ``+`` ``A`` ``-`` `[`min`](https://rdrr.io/r/base/Extremes.html)`(``Y``)``)`` ``/`` `[`diff`](https://rdrr.io/r/base/diff.html)`(`[`range`](https://rdrr.io/r/base/range.html)`(``Y``)``)`\
`Qn_upshift`` ``<-`` ``(``W`` ``+`` ``A`` ``+`` ``delta`` ``-`` `[`min`](https://rdrr.io/r/base/Extremes.html)`(``Y``)``)`` ``/`` `[`diff`](https://rdrr.io/r/base/diff.html)`(`[`range`](https://rdrr.io/r/base/range.html)`(``Y``)``)`\
`Qn_noshift``[``Qn_noshift`` ``<`` ``0``]`` ``<-`` ``0.025`\
`Qn_noshift``[``Qn_noshift`` ``>`` ``1``]`` ``<-`` ``0.975`\
`Qn_upshift``[``Qn_upshift`` ``<`` ``0``]`` ``<-`` ``0.025`\
`Qn_upshift``[``Qn_upshift`` ``>`` ``1``]`` ``<-`` ``0.975`\
`Qn_out`` ``<-`` `[`as.data.table`](https://rdrr.io/pkg/data.table/man/as.data.table.html)`(`[`cbind`](https://rdrr.io/r/base/cbind.html)`(``Qn_noshift``, ``Qn_upshift``)``)``[``C`` ``==`` ``1``, ``]`\
[`setnames`](https://rdrr.io/pkg/data.table/man/setattr.html)`(``Qn_out``, `[`c`](https://rdrr.io/r/base/c.html)`(``"noshift"``, ``"upshift"``)``)`\
\
`# construct efficient estimator by applying wrapper function`\
`est_shift_spec`` ``<-`` `[`txshift`](https://codex.nimahejazi.org/txshift/reference/txshift.md)`(`\
`  W ``=`` ``W``, A ``=`` ``A``, Y ``=`` ``Y``, delta ``=`` ``delta``,`\
`  samp_fit_args ``=`` ``NULL``,`\
`  samp_fit_ext ``=`` ``ipcw_out``,`\
`  g_exp_fit_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``fit_type ``=`` ``"external"``)``,`\
`  Q_fit_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``fit_type ``=`` ``"external"``)``,`\
`  gn_exp_fit_ext ``=`` ``gn_out``,`\
`  Qn_fit_ext ``=`` ``Qn_out`\
`)`

## References

Coyle, Jeremy R, Nima S Hejazi, Ivana Malenica, Rachael V Phillips, and
Oleg Sofrygin. 2022. *sl3: Modern Machine Learning Pipelines for Super
Learning*. Released. <https://doi.org/10.5281/zenodo.1342293>.

Coyle, Jeremy R, Nima S Hejazi, Rachael V Phillips, Lars W van der Laan,
and Mark J van der Laan. 2022. *hal9001: The Scalable Highly Adaptive
Lasso*. Released. <https://doi.org/10.5281/zenodo.3558313>.

Dı́az, Iván, and Mark J van der Laan. 2012. “Population Intervention
Causal Effects Based on Stochastic Interventions.” *Biometrics* 68 (2):
541–49.

Dı́az, Iván, and Mark J van der Laan. 2018. “Stochastic Treatment
Regimes.” In *Targeted Learning in Data Science: Causal Inference for
Complex Longitudinal Studies*. Springer Science & Business Media.

Dudoit, Sandrine, and Mark J van der Laan. 2005. “Asymptotics of
Cross-Validated Risk Estimation in Estimator Selection and Performance
Assessment.” *Statistical Methodology* 2 (2): 131–54.

Gilbert, Peter B, Youyi Fong, David Benkeser, et al. 2021. *CoVPN
COVID-19 Vaccine Efficacy Trial Immune Correlates Statistical Analysis
Plan*. <https://doi.org/10.6084/m9.figshare.13198595>.

Haneuse, Sebastian, and Andrea Rotnitzky. 2013. “Estimation of the
Effect of Interventions That Modify the Received Treatment.” *Statistics
in Medicine* 32 (30): 5260–77.

Hejazi, Nima S, Mark J van der Laan, Holly E Janes, Peter B Gilbert, and
David C Benkeser. 2020. “Efficient Nonparametric Inference on the
Effects of Stochastic Interventions Under Two-Phase Sampling, with
Applications to Vaccine Efficacy Trials.” *Biometrics* 77 (4): 1241–53.
<https://doi.org/10.1111/biom.13375>.

Pearl, Judea. 2000. *Causality*. Cambridge University Press.

Rose, Sherri, and Mark J van der Laan. 2011. “A Targeted Maximum
Likelihood Estimator for Two-Stage Designs.” *International Journal of
Biostatistics* 7 (1): 1–21.

van der Laan, Mark J, Jeremy R Coyle, Nima S Hejazi, Ivana Malenica,
Rachael V Phillips, and Alan E Hubbard. 2022. *Targeted Learning in `R`:
A Causal Data Science Handbook*. CRC Press.
<https://tlverse.org/tlverse-handbook/>.

van der Laan, Mark J, Sandrine Dudoit, and Sunduz Keles. 2004.
“Asymptotic Optimality of Likelihood-Based Cross-Validation.”
*Statistical Applications in Genetics and Molecular Biology* 3 (1):
1–23.

van der Laan, Mark J, Eric C Polley, and Alan E Hubbard. 2007. “Super
Learner.” *Statistical Applications in Genetics and Molecular Biology* 6
(1).

van der Laan, Mark J, and Sherri Rose. 2011. *Targeted Learning: Causal
Inference for Observational and Experimental Data*. Springer Science &
Business Media.

van der Laan, Mark J, and Sherri Rose. 2018. *Targeted Learning in Data
Science: Causal Inference for Complex Longitudinal Studies*. Springer
Science & Business Media.

van der Vaart, Aad W, Sandrine Dudoit, and Mark J van der Laan. 2006.
“Oracle Inequalities for Multi-Fold Cross Validation.” *Statistics &
Decisions* 24 (3): 351–71.
