#' Method to get the sojourn time distribution f
#' 
#' @description Computes the conditional sojourn time distribution \eqn{f(k)}, 
#'   \eqn{f_{i}(k)}, \eqn{f_{j}(k)} or \eqn{f_{ij}(k)}.
#'   
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric] or [smmnonparametric]).
#' @param k A positive integer giving the time horizon.
#' @return A vector, matrix or array giving the value of \eqn{f} at each time 
#'   between 0 and `k`.
#'   
#' @noRd
#' 
.get.f <- function(x, k) {
  UseMethod(".get.f", x)
}


#' Method to get the number of parameters of the semi-Markov chain
#' 
#' @description Method to get the number of parameters of the semi-Markov 
#'   chain. This method is useful for the computation of criteria such as AIC 
#'   and BIC.
#' 
#' @param x An object for which the number of parameters can be returned.
#' @return A positive integer giving the number of parameters.
#' 
#' @noRd
#' 
.getKpar <- function(x) {
  UseMethod(".getKpar", x)
}


#' Log-likelihood Function
#' 
#' @description Computation of the log-likelihood for a semi-Markov model
#' 
#' @param x An object for which the log-likelihood can be computed.
#' @param processes An object of class `processes`.
#' 
#' @noRd
#' 
.loglik <- function(x, processes) {
  UseMethod(".loglik", x)
}


#' Akaike Information Criterion (AIC)
#' 
#' @description Generic function computing the Akaike Information Criterion of 
#'   the model `x`, with the list of sequences `sequences`.
#' 
#' @param x An object for which there exists a `loglik` method to compute the 
#'   log-likelihood.
#' @param sequences A list of vectors representing the sequences for which the 
#'   AIC will be computed based on `x`.
#' @return Value of the AIC.
#' 
#' @export
#' 
aic <- function(x, sequences) {
  UseMethod("aic", x)
}


#' Bayesian Information Criterion (BIC)
#' 
#' @description Generic function computing the Bayesian Information Criterion 
#'   of the model `x`, with the list of sequences `sequences`.
#' 
#' @param x An object for which there exists a `loglik` method to compute the 
#'   log-likelihood.
#' @param sequences A list of vectors representing the sequences for which the 
#'   BIC will be computed based on `x`.
#' @return Value of the BIC.
#' 
#' @export
#' 
bic <- function(x, sequences) {
  UseMethod("bic", x)
}


#' Method to get the semi-Markov kernel \eqn{q}
#' 
#' @description Computes the semi-Markov kernel \eqn{q_{ij}(k)}.
#' 
#' @param x An object inheriting from the S3 class `smm` (an object of class
#'   [smmparametric] or [smmnonparametric]).
#' @param k A positive integer giving the time horizon.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes] function) for the asymptotic variance.
#' @return An array giving the value of \eqn{q_{ij}(k)} at each time between 0 
#'   and `k` if `var = FALSE`. If `var = TRUE`, a list containing the 
#'   following components:
#'   \itemize{
#'    \item{x: }{an array giving the value of \eqn{q_{ij}(k)} at each time 
#'      between 0 and `k`;}
#'    \item{sigma2: }{an array giving the asymptotic variance of the estimator 
#'      \eqn{\sigma_{q}^{2}(i, j, k)}.}
#'  }
#'  
#' @export
#' 
getKernel <- function(x, k, var = FALSE, klim = 10000) {
  UseMethod("getKernel", x)
}


#' Log-likelihood Function
#' 
#' @description Generic function computing the log-likelihood of the model `x`,
#'   with the list of sequences `sequences`.
#' 
#' @param x An object for which the log-likelihood can be computed.
#' @param sequences A list of vectors representing the sequences for which the 
#'   log-likelihood will be computed based on `x`.
#' @return Value of the log-likelihood.
#' 
#' @export
#' 
loglik <- function(x, sequences) {
  UseMethod("loglik", x)
}
