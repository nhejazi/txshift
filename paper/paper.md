---
title: "`txshift`: Efficient estimation of the causal effects of stochastic interventions in `R`"
tags:
  - causal inference
  - machine learning
  - two-phase sampling
  - efficient estimation
  - targeted learning
  - stochastic intervention
  - modified treatment policy
  - R
authors:
  - name: Nima S. Hejazi
    orcid: 0000-0002-7127-2789
    affiliation: 1, 2
  - name: David Benkeser
    orcid: 0000-0002-1019-8343
    affiliation: 3
affiliations:
  - name: Graduate Group in Biostatistics, University of California, Berkeley
    index: 1
  - name: Center for Computational Biology, University of California, Berkeley
    index: 2
  - name: Department of Biostatistics and Bioinformatics, Rollins School of Public Health, Emory University
    index: 3
date: 28 September 2020
bibliography: refs.bib
---

# Summary

The `txshift` `R` package aims to provide researchers in (bio)statistics,
epidemiology, health policy, econometrics, and related disciplines with access
to state-of-the-art statistical methodology for evaluating the causal effects of
stochastic shift interventions on _continuous-valued_ exposures. `txshift`
estimates the causal effects of modified treatment policies (or "feasible
interventions"), which take into account the natural value of an exposure in
assigning an intervention level. To accommodate use in study designs
incorporating two-phase sampling (e.g., case-control), the package provides two
types of modern corrections, both rooted in semiparametric theory, for
constructing unbiased and efficient estimates, despite the significant
limitations induced by such designs. Thus, `txshift` makes possible the
estimation of the causal effects of stochastic interventions in experimental and
observational study settings subject to real-world design limitations that
commonly arise in modern scientific practice.

# Background

Causal inference has traditionally focused on the effects of static
interventions, under which the magnitude of the exposure is set to a fixed,
prespecified value for each unit. The evaluation of such interventions faces
a host of issues, among them non-identification, violations of the assumption of
positivity, and inefficiency. Stochastic interventions provide a promising
solution to these fundamental issues by allowing for the target parameter to be
defined as the mean counterfactual outcome under a hypothetically shifted
version of the observed exposure distribution [@diaz2012population].
Modified treatment policies, a particular class of such interventions, may be
interpreted as shifting the natural exposure level at the level of a given
observational unit [@haneuse2013estimation; @diaz2018stochastic].

Despite the promise of such advances in causal inference, real data analyses are
often further complicated by economic constraints, such as when the primary
variable of interest is far more expensive to collect than auxiliary covariates.
Two-phase sampling is often used to bypass these limitations -- unfortunately,
these sampling schemes produce side effects that require further adjustment when
formal statistical inference is the principal goal of a study. Among the rich
literature on two-phase designs, @rose2011targeted2sd stand out for providing
a study of nonparametric efficiency theory under a broad class of two-phase
designs. Their work provides guidance on constructing efficient estimators of
causal effects under general two-phase sampling designs.

# `txshift`'s Scope

Building on these prior works, @hejazi2020efficient outlined a novel approach
for use in such settings: augmented targeted minimum loss (TML) and one-step
estimators for the causal effects of stochastic interventions, with guarantees
of consistency, efficiency, and multiple robustness even in the presence of
two-phase sampling. These authors further outlined a technique that summarizes
the effect of shifting an exposure variable on the outcome of interest via
a nonparametric working marginal structural model, analogous to a dose-response
analysis. The `txshift` software package, for the `R` language and environment
for statistical computing [@R], implements this methodology.

`txshift` is designed to facilitate the construction of TML and one-step
estimators of the causal effects of modified treatment policies that shift the
observed exposure value up (or down) by an arbitrary scalar $\delta$, which may
possibly take into account the natural value of the exposure (and, in future
versions, the covariates). The `R` package includes tools for deploying these
efficient estimators under two-phase sampling designs, with two types of
corrections: (1) a reweighting procedure that introduces inverse probability of
censoring weights directly into relevant loss functions, as discussed in
@rose2011targeted2sd; as well as (2) an augmented efficient influence function
estimating equation, studied more thoroughly by @hejazi2020efficient. `txshift`
integrates with the [`sl3` package](https://github.com/tlverse/sl3)
[@coyle2020sl3] to allow for ensemble machine learning to be leveraged in the
estimation of nuisance parameters. What's more, the `txshift` package draws on
both the `hal9001` [@coyle2020hal9001; @hejazi2020hal9001] and `haldensify`
[@hejazi2020haldensify] `R` packages to allow each of the efficient estimators
to be constructed in a manner consistent with the methodological and theoretical
advances of @hejazi2020efficient, which require fast convergence rates of
nuisance parameters to their true counterparts for efficiency of the resultant
estimator.

# Availability

The `txshift` package has been made publicly available both [via
GitHub](https://github.com/nhejazi/txshift) and the [Comprehensive `R` Archive
Network](https://CRAN.R-project.org/package=txshift). Use of the `txshift`
package has been extensively documented in the package's `README`, two
vignettes, and its [`pkgdown` documentation
website](https://code.nimahejazi.org/txshift).

# Acknowledgments

Nima Hejazi's contributions to this work were supported in part by a grant from
the National Institutes of Health: [T32
LM012417-02](https://projectreporter.nih.gov/project_info_description.cfm?aid=9248418&icde=37849831&ddparam=&ddvalue=&ddsub=&cr=1&csb=default&cs=ASC&pball=).

# References

