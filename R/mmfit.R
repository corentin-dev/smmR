mmfit <- function(mm, M, logLik, sequences) {
  
  ans <- mm
  ans$M <- M
  ans$logLik <- logLik
  ans$sequences <- sequences
  
  class(ans) <- c("mmfit", class(ans))
  
  return(ans)
  
}


#' Function to check if an object is of class `mmfit`
#' 
#' @description `is.mmfit` returns `TRUE` if `x` is an object of 
#'   class `mmfit`.
#' 
#' @param x An arbitrary R object.
#' @return `is.mmfit` returns `TRUE` or `FALSE` depending on whether `x` is an 
#'   object of class `mmfit` or not.
#' 
#' @export
#' 
is.mmfit <- function(x) {
  inherits(x, "mmfit")
}


# Method to get the number of parameters
# (useful for the computation of criteria such as AIC and BIC)
.get.Kpar.mmfit <- function(x) {
  NextMethod(x)
}


#' Akaike Information Criterion (AIC)
#' 
#' @description Computation of the Akaike Information Criterion.
#' 
#' @param x An object of class `mmfit`.
#' @param sequences A list of vectors representing the sequences for which the 
#'   AIC will be computed based on `x`.
#' @return Value of the AIC.
#' 
#' @noRd
#' 
#' @export
#' 
AIC.mmfit <- function(object, ...) {
  
  sequences = list(...)[1]
  logLik <- logLik(object, sequences)
  
  kpar <- .get.Kpar(object)
  
  AIC <- -2 * logLik + 2 * kpar
  
  return(AIC)
  
}


#' Bayesian Information Criterion (BIC)
#' 
#' @description Computation of the Bayesian Information Criterion.
#' 
#' @param x An object of class `mmfit`.
#' @param sequences A list of vectors representing the sequences for which the 
#'   BIC will be computed based on `x`.
#' @return Value of the BIC.
#' 
#' @noRd
#' 
#' @export
#' 
BIC.mmfit <- function(object, ...) {
  
  sequences = list(...)[1]
  logLik <- logLik(object, sequences)
  
  kpar <- .get.Kpar(object)
  
  if (is.null(sequences)) {
    n <- object$M
  } else {
    n <- sum(sapply(sequences, length))  
  }
  
  BIC <- -2 * logLik + log(n) * kpar
  
  return(BIC)
  
}


#' Log-likelihood Function
#' 
#' @description Computation of the log-likelihood for a Markov model.
#' 
#' @param x An object of class `mmfit`.
#' @param sequences A list of vectors representing the sequences for which the 
#'   log-likelihood will be computed based on `x`.
#' @return Value of the log-likelihood.
#' 
#' @noRd
#' 
#' @export
#' 
logLik.mmfit <- function(object, ...) {
  
  # Computing a new value of log-likelihood based on the parameter sequences
  sequences = list(...)[1]
  if (!is.null(sequences)) {
    
    if (!(is.list(sequences) & all(sapply(sequences, class) %in% c("character", "numeric")))) {
      stop("The parameter 'sequences' should be a list of vectors")
    }
    
    if (!all(unique(unlist(sequences)) %in% object$states)) {
      stop("Some states in the list of observed sequences 'sequences' 
         are not in the state space given by the model 'object'")
    }
    
    NextMethod(object)
    
  } else {# Return the value of the log-likelihood
    
    return(object$logLik)
    
  }
  
}


#' Simulates Markov chains
#' 
#' @description Simulates sequences from a fitted Markov model.
#' 
#' @details If `nsim` is a single integer then a chain of that length is 
#'   produced. If `nsim` is a vector of integers, then `length(nsim)` 
#'   sequences are generated with respective lengths.
#' 
#' @param object An object of class `mmfit`.
#' @param nsim An integer or vector of integers (for multiple sequences) 
#'   specifying the length of the sequence(s).
#' @param seed Optional. `seed` for the random number generator. 
#'   If no `seed` is given, then seed is set by using the command 
#'   `set.seed(round(as.numeric(Sys.time()))`.
#' @param ... further arguments passed to or from other methods.
#' @return A list of vectors representing the sequences.
#' 
#' @seealso [mm], [fitmm]
#' 
#' @export
#' 
simulate.mmfit <- function(object, nsim = 1, seed = NULL, ...) {
  NextMethod(object)
}
