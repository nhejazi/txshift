---
title: "`txshift`: Efficient estimation of the causal effects of stochastic interventions in `R`"
tags:
  - causal inference
  - machine learning
  - two-phase sampling
  - targeted learning
  - stochastic interventions
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
date: 30 March 2020
bibliography: refs.bib
---

# Summary

Causal inference has traditionally focused on the effects of static
interventions, under which the magnitude of the treatment is set to a fixed,
prespecified value for each unit. The evaluation of such interventions faces
a host of issues, among them non-identification, violations of the assumption of
positivity, and inefficiency. Stochastic interventions provide a promising
solution to these fundamental issues by allowing for the target parameter to be
defined as the mean counterfactual outcome under a hypothetically shifted
version of the observed exposure distribution [@diaz2012population].
@rotnitzky2013... developed...
Subsequent
work by @diaz2018stochastic revealed a simplified estimation strategy for
estimating the effects

Despite the promise of such approaches, real data analyses are often further
complicated by economic constraints, such as when the primary variable of
interest is far more expensive to collect than auxiliary covariates. Two-phase
sampling schemes are often used to bypass such limitations -- unfortunately,
their use produces side effects that require further adjustment when formal
statistical inference is the principal goal of a study. While a rich history of
work on estimation and inference under such designs has developed,
@rose2011targeted2sd provided a unique treatment of the subject, formulating...

Building on these prior works, @hejazi2020efficient outlined a novel approach
for use in such settings: augmented targeted minimum loss (TML) and one-step
estimators for the causal effects of stochastic interventions, with guarantees
of consistency, efficiency, and multiple robustness even in the presence of
two-phase sampling. These authors further proposed a technique utilizing the
estimated causal effects of stochastic interventions to construct
a nonparametric working marginal structural model to summarize the effect of
shifting an exposure variable on the outcome of interest, analogous to
a dose-response analysis. The `txshift` software package, for the `R` language
and environment for statistical computing [@R], is an implementation of this
methodology.

`txshift` provides...


# Acknowledgments

Nima Hejazi's contributions to this work were supported in part by a grant from
the National Institutes of Health: [T32
LM012417-02](https://projectreporter.nih.gov/project_info_description.cfm?aid=9248418&icde=37849831&ddparam=&ddvalue=&ddsub=&cr=1&csb=default&cs=ASC&pball=).

# References

