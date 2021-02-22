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
#' Semi-Markov models are specified by using the functions `smmparametric()` 
#' and `smmnonparametric()` for parametric and non-parametric specifications 
#' respectively. These functions return objects of S3 class (`smm`, 
#' `smmparametric`) and (`smm`, `smmnonparametric`) respectively (`smm` class 
#' inherits from S3 classes `smmparametric` or `smmnonparametric`). Thus, `smm` 
#' is like a wrapper class for semi-Markov model specifications.
#' 
#' Based on a model specification (an object of class `smm`), it is possible to:
#' \itemize{
#'   \item **simulate** one or several sequences with the method 
#'     `simulate.smm()`;
#'   \item **plot** conditional sojourn time distributions (method 
#'     `plot.smm()`);
#'   \item compute **log-likelihood**, **AIC** and **BIC** criteria (methods 
#'     `loglik()`, `aic()`, `bic()`);
#'   \item compute **reliability**, **maintainability**, **availability**, 
#'     **failure rates** (methods `reliability()`, `maintainability()`, 
#'     `availability()`, `failureRate()`).
#' }
#' 
#' Estimations of parametric and non-parametric semi-Markov models can be done 
#' by using the function `fitsemimarkovmodel()`. This function returns an 
#' object of S3 class `fitsmm`. `fitsmm` inherits from classes 
#' (`smm`, `smmparametric`) or (`smm`, `smmnonparametric`).
#' 
#' Based on a fitted/estimated semi-Markov model (an object of class `smmfit`),
#' it is possible to:
#' \itemize{
#'   \item **simulate** one or several sequences with the method 
#'     `simulate.smmfit()`;
#'   \item **plot** estimated conditional sojourn time distributions 
#'     (method `plot.smmfit()`);
#'   \item compute **log-likelihood**, **AIC** and **BIC** criteria (methods 
#'     `loglik()`, `aic()`, `bic()`);
#'   \item compute estimated **reliability**, **maintainability**, 
#'     **availability**, **failure rates** and their **confidence intervals** 
#'     (methods `reliability()`, `maintainability()`, `availability()`, 
#'     `failureRate()`).
#' }
#' 
#' @keywords markov semi-markov simulation estimation censored
#' 
#' @import utils
#' @import stats
#' @import seqinr
#' @import DiscreteWeibull
#' 
## usethis namespace: start
#' @useDynLib SMM, .registration = TRUE
## usethis namespace: end
## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
"_PACKAGE"
