smmfit <- function(smm, M, loglik, sequences) {
  
  ans <- smm
  ans$M <- M
  ans$loglik <- loglik
  ans$sequences <- sequences
  
  class(ans) <- c("smmfit", class(ans))
  
  return(ans)
  
}


#' Function to check if an object is of class `smmfit`
#' 
#' @description `is.smmfit` returns `TRUE` if `x` is an object of class `smmfit`.
#' 
#' @param x An arbitrary R object.
#' @return `is.smmfit` returns `TRUE` or `FALSE` depending on whether `x` is an 
#'   object of class `smmfit` or not.
#' 
#' @export
#' 
is.smmfit <- function(x) {
  inherits(x, "smmfit")
}


# Method to get the number of parameters
# (useful for the computation of criteria such as AIC and BIC)
.getKpar.smmfit <- function(x) {
  NextMethod(x)
}


#' Akaike Information Criterion (AIC)
#' 
#' @description Computation of the Akaike Information Criterion.
#' 
#' @param x An object of class `smmfit` (inheriting from the S3 classes 
#'   `smm`, [smmnonparametric] or [smmparametric]).
#' @param sequences A list of vectors representing the sequences for which the 
#'   AIC will be computed based on `x`.
#' @return Value of the AIC.
#' 
#' @noRd
#' 
#' @export
#' 
aic.smmfit <- function(x, sequences = NULL) {
  
  loglik <- loglik(x, sequences)
  
  kpar <- .getKpar(x)
  
  aic <- -2 * loglik + 2 * kpar
  
  return(aic)
  
}


#' Bayesian Information Criterion (BIC)
#' 
#' @description Computation of the Bayesian Information Criterion.
#' 
#' @param x An object of class `smmfit` (inheriting from the S3 classes 
#'   `smm`, [smmnonparametric] or [smmparametric]).
#' @param sequences A list of vectors representing the sequences for which the 
#'   BIC will be computed based on `x`.
#' @return Value of the BIC.
#' 
#' @noRd
#' 
#' @export
#' 
bic.smmfit <- function(x, sequences = NULL) {
  
  loglik <- loglik(x, sequences)
  
  kpar <- .getKpar(x)
  
  if (is.null(sequences)) {
    n <- x$M
  } else {
    n <- sum(sapply(sequences, length))  
  }
  
  bic <- -2 * loglik + log(n) * kpar
  
  return(bic)
  
}


#' Method to get the semi-Markov kernel \eqn{q}
#' 
#' @description Computes the semi-Markov kernel \eqn{q_{ij}(k)}.
#' 
#' @param x An object of class `smmfit` (inheriting from the S3 classes 
#'   `smm`, [smmnonparametric] or [smmparametric]).
#' @param k A positive integer giving the time horizon.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param klim Optional. The time horizon used to approximate the series in the
#'   computation of the mean recurrence times vector for the asymptotic 
#'   variance.
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
#' @noRd
#' 
#' @export
#' 
getKernel.smmfit <- function(x, k, var = FALSE, klim = 10000) {
  NextMethod(x)
}


#' Log-likelihood Function
#' 
#' @description Computation of the log-likelihood for a semi-Markov model.
#' 
#' @param x An object of class `smmfit` (inheriting from the S3 classes 
#'   `smm`, [smmnonparametric] or [smmparametric]).
#' @param sequences A list of vectors representing the sequences for which the 
#'   log-likelihood will be computed based on `x`.
#' @return Value of the log-likelihood.
#' 
#' @noRd
#' 
#' @export
#' 
loglik.smmfit <- function(x, sequences = NULL) {
  
  # Computing a new value of log-likelihood based on the parameter sequences
  if (!is.null(sequences)) {
    
    if (!(is.list(sequences) & all(sapply(sequences, class) %in% c("character", "numeric")))) {
      stop("The parameter 'sequences' should be a list of vectors")
    }
    
    if (!all(unique(unlist(sequences)) %in% x$states)) {
      stop("Some states in the list of observed sequences 'sequences' 
         are not in the state space given by the model 'x'")
    }
    
    NextMethod(x)
    
  } else {# Return the value of the log-likelihood
    
    return(x$loglik)
    
  }
  
}


#' Plot function for an object of class smmfit
#' 
#' @description Displays the densities for the conditional sojourn time 
#'   distributions depending on the current state `i` and on the next state 
#'   `j`.
#'   
#' @param x An object of class `smmfit` (inheriting from the S3 classes 
#'   `smm`, [smmnonparametric] or [smmparametric]).
#' @param i An element of the state space vector `x$states` giving the current 
#'   state in the following cases: `type.sojourn = "fij"` or `type.sojourn = "fi"`, 
#'   otherwise, `i` is ignored.
#' @param j An element of the state space vector `x$states` giving the next 
#'   state in the following cases: `type.sojourn = "fij"` or `type.sojourn = "fj"`, 
#'   otherwise, `j` is ignored.
#' @param klim An integer giving the limit value for which the density will be 
#'   plotted. If `klim` is `NULL`, then quantile or order 0.95 is used.
#' @param ... Arguments passed to plot.
#' @return None.
#' 
#' @references
#' V. S. Barbu, N. Limnios. (2008). Semi-Markov Chains and Hidden Semi-Markov 
#' Models Toward Applications - Their Use in Reliability and DNA Analysis. 
#' New York: Lecture Notes in Statistics, vol. 191, Springer.
#' 
#' @export
#' 
#' @examples 
#' states <- c("a", "c", "g", "t")
# s <- length(states)
# 
# # Creation of the initial distribution
# vect.init <- c(1 / 4, 1 / 4, 1 / 4, 1 / 4)
# 
# # Creation of the transition matrix
# pij <- matrix(c(0, 0.2, 0.3, 0.4, 0.2, 0, 0.5, 0.2, 0.5, 
#                 0.3, 0, 0.4, 0.3, 0.5, 0.2, 0), ncol = s)
# 
# # Creation of the distribution matrix
# distr.vec <- c("pois", "geom", "geom", "geom")
# parameters <- matrix(c(2, 0.6, 0.8, 0.8, NA, NA, NA, NA), 
#                      ncol = 2, byrow = FALSE)
# 
# # Specify the semi-Markov model
# smm <- smmparametric(states = states, init = vect.init, ptrans = pij, 
#                      type.sojourn = "fi", distr = distr.vec, param = parameters)
# 
# seqs <- simulate(object = smm, nsim = rep(5000, 100), seed = 10)
# 
# est <- fitsmm(sequences = seqs, states = states, type.sojourn = "fi", distr = distr.vec)
# 
# class(est)
# 
# plot(x = est, i = "a", col = "blue", pch = "+")
# 
plot.smmfit <- function(x, i, j, klim = NULL, ...) {
  NextMethod(x)
}


#' Simulates semi-Markov chains
#' 
#' @description Simulates sequences from a fitted semi-Markov model.
#' 
#' @details If `nsim` is a single integer then a chain of that length is 
#'   produced. If `nsim` is a vector of integers, then `length(nsim)` 
#'   sequences are generated with respective lengths.
#' 
#' @param object An object of class `smmfit` (inheriting from the S3 classes 
#'   `smm`, [smmnonparametric] or [smmparametric]).
#' @param nsim An integer or vector of integers (for multiple sequences) 
#'   specifying the length of the sequence(s).
#' @param seed `seed` for the random number generator.
#' @param ... further arguments passed to or from other methods.
#' @return A list of vectors representing the sequences.
#' 
#' @seealso [smmnonparametric], [smmparametric], [fitsmm]
#' 
#' @references
#' V. S. Barbu, N. Limnios. (2008). Semi-Markov Chains and Hidden Semi-Markov 
#' Models Toward Applications - Their Use in Reliability and DNA Analysis. 
#' New York: Lecture Notes in Statistics, vol. 191, Springer.
#' 
#' @export
#' 
simulate.smmfit <- function(object, nsim = 1, seed = NULL, ...) {
  NextMethod(object)
}


#' @export
reliability.smmfit <- function(x, k, upstates = x$states, level = 0.95, klim = 10000) {
  
  out <- NextMethod()
  
  out <- cbind(out, out[, 1] + (sqrt(out[, 2] / x$M) * qnorm(p = 1 - (1 - level) / 2)) %o% c(-1, 1))
  colnames(out) <- c("reliability", "sigma2", "lwr", "upper")
  
  return(out)
}


#' @export
maintainability.smmfit <- function(x, k, upstates = x$states, level = 0.95, klim = 10000) {
  
  out <- NextMethod()
  
  out <- cbind(out, out[, 1] + (sqrt(out[, 2] / x$M) * qnorm(p = 1 - (1 - level) / 2)) %o% c(-1, 1))
  colnames(out) <- c("maintainability", "sigma2", "lwr", "upper")
  
  return(out)
}


#' @export
availability.smmfit <- function(x, k, upstates = x$states, level = 0.95, klim = 10000) {
  
  out <- NextMethod()
  
  out <- cbind(out, out[, 1] + (sqrt(out[, 2] / x$M) * qnorm(p = 1 - (1 - level) / 2)) %o% c(-1, 1))
  colnames(out) <- c("availability", "sigma2", "lwr", "upper")
  
  return(out)
}


.failureRateBMP.smmfit <- function(x, k, upstates = x$states, level = 0.95, epsilon = 1e-3, klim = 10000) {
  
  out <- .failureRateBMP.smm(x = x, k = k, upstates = upstates, level = level, epsilon = epsilon, klim = klim)
  
  out <- cbind(out, out[, 1] + (sqrt(out[, 2] / x$M) * qnorm(p = 1 - (1 - level) / 2)) %o% c(-1, 1))
  colnames(out) <- c("BMP", "sigma2", "lwr", "upper")
  
  return(out)
}


.failureRateRG.smmfit <- function(x, k, upstates = x$states, level = 0.95, epsilon = 1e-3, klim = 10000) {
  
  out <- .failureRateRG.smm(x = x, k = k, upstates = upstates, level = level, epsilon = epsilon, klim = klim)
  
  out <- cbind(out, out[, 1] + (sqrt(out[, 2] / x$M) * qnorm(p = 1 - (1 - level) / 2)) %o% c(-1, 1))
  colnames(out) <- c("RG", "sigma2", "lwr", "upper")
  
  return(out)
}


#' @export
failureRate.smmfit <- function(x, k, upstates = x$states, failure.rate = c("BMP", "RG"), level = 0.95, epsilon = 1e-3, klim = 10000) {
  
  #############################
  # Checking parameters failure.rate
  #############################
  
  failure.rate <- match.arg(failure.rate)
  
  
  if (failure.rate == "BMP") {
    out <- .failureRateBMP.smmfit(x = x, k = k, upstates = upstates, level = level, epsilon = epsilon, klim = klim)
  } else {
    out <- .failureRateRG.smmfit(x = x, k = k, upstates = upstates, level = level, epsilon = epsilon, klim = klim)
  }
  
  return(out)
  
}


#' @export
mttf.smmfit <- function(x, upstates = x$states, level = 0.95, klim = 10000) {
  
  out <- NextMethod()
  
  out <- cbind(out, out[, 1] + (sqrt(out[, 2] / x$M) * qnorm(p = 1 - (1 - level) / 2)) %o% c(-1, 1))
  colnames(out) <- c("mttf", "sigma2", "lwr", "upper")
  
  return(out)
}


#' @export
mttr.smmfit <- function(x, upstates = x$states, level = 0.95, klim = 10000) {
  
  out <- NextMethod()
  
  out <- cbind(out, out[, 1] + (sqrt(out[, 2] / x$M) * qnorm(p = 1 - (1 - level) / 2)) %o% c(-1, 1))
  colnames(out) <- c("mttr", "sigma2", "lwr", "upper")
  
  return(out)
}
