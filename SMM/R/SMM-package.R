#' SMM : Semi-Markov and Markov Models
#'
#' @description This package performs parametric and non-parametric estimation 
#' and simulation for multi-state discrete-time semi-Markov processes. For the 
#' parametric estimation, several discrete distributions are considered for the
#' sojourn times: Uniform, Geometric, Poisson, Discrete Weibull and Negative 
#' Binomial. The non-parametric estimation concerns the sojourn time 
#' distributions, where no assumptions are done on the shape of distributions. 
#' Moreover, the estimation can be done on the basis of one or several sample 
#' paths, with or without censoring at the beginning or/and at the end of the 
#' sample paths. Estimation and simulation of discrete-time k-th order Markov 
#' chains are also considered.
#'
#' @keywords markov semi-markov simulation estimation censored
#'
#' @import utils
#' @import stats
#' @import seqinr
#' @import DiscreteWeibull
"_PACKAGE"