markovmodelfit <- function(markovmodel, M, loglik, sequences) {
  
  ans <- markovmodel
  ans$M <- M
  ans$loglik <- loglik
  ans$sequences <- sequences
  
  class(ans) <- c("markovmodelfit", class(ans))
  
  return(ans)
  
}


#' Function to check if an object is of class `markovmodelfit`
#' 
#' @description `is.markovmodelfit` returns `TRUE` if `x` is an object of 
#'   class `markovmodelfit`.
#' 
#' @param x An arbitrary R object.
#' 
#' @export
#' 
is.markovmodelfit <- function(x) {
  inherits(x, "markovmodelfit")
}


# Method to get the number of parameters
# (useful for the computation of criteria such as AIC and BIC)
.getKpar.markovmodelfit <- function(x) {
  NextMethod(x)
}


#' Akaike Information Criterion (AIC)
#' 
#' @description Computation of the Akaike Information Criterion.
#' 
#' @param x An object of class `markovmodelfit`.
#' @param sequences A list of vectors representing the sequences for which the 
#'   AIC will be computed based on `x`.
#' @return Value of the AIC.
#' 
#' @noRd
#' 
#' @export
#' 
aic.markovmodelfit <- function(x, sequences = NULL) {
  
  loglik <- loglik(x, sequences)
  
  kpar <- .getKpar(x)
  
  aic <- -2 * loglik + 2 * kpar
  
  return(aic)
  
}


#' Bayesian Information Criterion (BIC)
#' 
#' @description Computation of the Bayesian Information Criterion.
#' 
#' @param x An object of class `markovmodelfit`.
#' @param sequences A list of vectors representing the sequences for which the 
#'   BIC will be computed based on `x`.
#' @return Value of the BIC.
#' 
#' @noRd
#' 
#' @export
#' 
bic.markovmodelfit <- function(x, sequences = NULL) {
  
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


#' Log-likelihood Function
#' 
#' @description Computation of the log-likelihood for a Markov model.
#' 
#' @param x An object of class `markovmodelfit`.
#' @param sequences A list of vectors representing the sequences for which the 
#'   log-likelihood will be computed based on `x`.
#' @return Value of the log-likelihood.
#' 
#' @noRd
#' 
#' @export
#' 
loglik.markovmodelfit <- function(x, sequences = NULL) {
  
  # Computing a new value of log-likelihood based on the parameter sequences
  if (!is.null(sequences)) {
    
    if (!(is.list(sequences) && all(sapply(sequences, class) %in% c("character", "numeric")))) {
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


#' Simulates Markov chains
#' 
#' @description Simulates sequences from a fitted Markov model.
#' 
#' @details If `nsim` is a single integer then a chain of that length is 
#'   produced. If `nsim` is a vector of integers, then `length(nsim)` 
#'   sequences are generated with respective lengths.
#' 
#' @param object An object of class `markovmodelfit`.
#' @param nsim An integer or vector of integers (for multiple sequences) 
#'   specifying the length of the sequence(s).
#' @param seed `seed` for the random number generator.
#' @param ... further arguments passed to or from other methods.
#' @return A list of vectors representing the sequences.
#' 
#' @seealso [markovmodel], [fitmarkovmodel]
#' 
#' @export
#' 
simulate.markovmodelfit <- function(object, nsim = 1, seed = NULL, ...) {
  NextMethod(object)
}
