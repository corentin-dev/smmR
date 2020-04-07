\name{SMM-package}
\alias{SMM-package}
\alias{SMM}
\docType{package}
\title{
SMM : Semi-Markov and Markov Models
}
\description{
This package performs parametric and non-parametric estimation and simulation for multi-state discrete-time semi-Markov processes. For the parametric estimation, several discrete distributions are considered for the sojourn times: Uniform, Geometric, Poisson, Discrete Weibull and Negative Binomial. The non-parametric estimation concerns the sojourn time distributions, where no assumptions are done on the shape of distributions. Moreover, the estimation can be done on the basis of one or several sample paths, with or without censoring at the beginning or/and at the end of the sample paths. Estimation and simulation of discrete-time k-th order Markov chains are also considered. 
}
\details{
This R package provides the following functions: estimSM, simulSM, LoglikelihoodSM, AIC_SM, BIC_SM, InitialLawSM, estimMk, simulMk, LoglikelihoodMK, AIC_Mk, BIC_Mk and InitialLawMk.
}

\author{
Vlad Stefan Barbu, barbu@univ-rouen.fr \cr Caroline Berard, caroline.berard@univ-rouen.fr \cr Dominique Cellier, dominique.cellier@laposte.net \cr Mathilde Sautreuil, mathilde.sautreuil@etu.univ-rouen.fr \cr Nicolas Vergne, nicolas.vergne@univ-rouen.fr 
}


\keyword{Markov models}
\keyword{Semi-Markov models }
\keyword{Simulation}
\keyword{Estimation}
\keyword{Censored data}
