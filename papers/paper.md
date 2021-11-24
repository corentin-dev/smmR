---
title: 'smmR: A Semi-Markov R package'
tags:
  - R
  - statistics
  - markov-chains
  - sequence
  - estimation
authors:
  - name: Vlad Stefan Barbu
    orcid: 0000-0002-0840-016X
    affiliation: 1
  - name: Caroline Berard
    orcid: 0000-0003-1386-488X
    affiliation: 1
  - name: Dominique Cellier
    affiliation: 1
  - name: Florian Lecocq
    affiliation: 1
  - name: Mathilde Sautreuil
    orcid: 0000-0003-1985-0307
    affiliation: 1
  - name: Corentin Lothode
    orcid: 0000-0002-8209-317X
    affiliation: 1
  - name: Nicolas Vergne
    affiliation: 1
affiliations:
 - name: University of Rouen Normandy
   index: 1
date: 22 November 2021
bibliography: paper.bib
---

# Summary

This package performs parametric and non-parametric estimation and simulation for multi-state discrete-time semi-Markov processes. For the parametric estimation, several discrete distributions are considered for the sojourn times: Uniform, Geometric, Poisson, Discrete Weibull and Negative Binomial. The non-parametric estimation concerns the sojourn time distributions, where no assumptions are done on the shape of distributions. Moreover, the estimation can be done on the basis of one or several sample paths, with or without censoring at the beginning or/and at the end of the sample paths. Estimation and simulation of discrete-time k-th order Markov chains are also considered.

Semi-Markov models are specified by using the functions `smmparametric()` and `smmnonparametric()` for parametric and non-parametric specifications respectively. These functions return objects of S3 class (`smm`, `smmparametric`) and (`smm`, `smmnonparametric`) respectively (`smm` class inherits from S3 classes `smmparametric` or `smmnonparametric`). Thus, `smm` is like a wrapper class for semi-Markov model specifications.

Based on a model specification (an object of class `smm`), it is possible to:

* simulate one or several sequences with the method `simulate.smm()`;
* plot conditional sojourn time distributions (method `plot.smm()`);
* compute log-likelihood, AIC and BIC criteria (methods `loglik()`, `aic()`, `bic()`);
* compute reliability, maintainability, availability, failure rates (methods `reliability()`, `maintainability()`, `availability()`, `failureRate()`).

Estimations of parametric and non-parametric semi-Markov models can be done by using the function `fitsmm()`. This function returns an object of S3 class `smmfit`. The class smmfit inherits from classes (`smm`, `smmparametric`) or (`smm`, `smmnonparametric`).

Based on a fitted/estimated semi-Markov model (an object of class `smmfit`), it is possible to:

* simulate one or several sequences with the method `simulate.smmfit()`;
* plot estimated conditional sojourn time distributions (method `plot.smmfit()`);
* compute log-likelihood, AIC and BIC criteria (methods `loglik()`, `aic()`, `bic()`);
* compute estimated reliability, maintainability, availability, failure rates and their confidence intervals (methods `reliability()`, `maintainability()`, `availability()`, `failureRate()`).

The implemented methods are described in:

* [@Barbu:2008]
* [@Barbu:2006]
* [@Trevezas:2011]

# Statement of need

Few R packages have been developed to handle semi-Markov models or hidden semi-Markov models. For semi-Markov models we have the recent semiMarkov R package [@Listwon:2015] that performs maximum likelihood estimation for parametric continuous-time semi-Markov processes, where the distribution can be chosen between Exponential, Weibull or exponentiated Weibull. That package computes associated hazard rates; covariates can also be taken into account through the Cox proportional hazard model. Two R packages are also dedicated to hidden semi- Markov models, implementing estimation and prediction methods: the hsmm R package [@Bulla:2010] and the mhsmm R package [@OConnell:2011]. Note that there is no R package developed for discrete-time multi-state semi-Markov models and none is computing reliability, maintainability, availability and failure rates. Thus the purpose of this package is to fill this need.

# Acknowledgements

We acknowledge financial contributions from

# References
