#' Markov model specification
#' 
#' @description Creates a model specification of a Markov model.
#' 
#' @param states Vector of state space of length s.
#' @param k Order of the Markov chain.
#' @param init Vector of initial distribution of length s ^ k.
#' @param ptrans Matrix of transition probabilities of dimension \eqn{(s, s)}.
#' @return An object of class [mm].
#' 
#' @seealso [simulate.mm], [fitmm]
#' 
#' @export
#' 
#' @examples
#' states <- c("a", "c", "g", "t")
#' s <- length(states)
#' k <- 1
#' init <- rep.int(1 / s, s)
#' p <- matrix(c(0, 0, 0.3, 0.4, 0, 0, 0.5, 0.2, 0.7, 0.5, 
#'               0, 0.4, 0.3, 0.5, 0.2, 0), ncol = s)
#' 
#' # Specify a Markov model of order 1
#' markov <- mm(states = states, init = init, ptrans = p, k = k)
#' 
mm <- function(states, init, ptrans, k = 1) {
  
  #############################
  # Checking parameter states
  #############################
  
  s <- length(states)
  
  if (!(is.vector(states) & (length(unique(states)) == s))) {
    stop("The state space 'states' is not a vector of unique elements")
  }
  
  #############################
  # Checking parameter init
  #############################
  
  if (!(is.numeric(init) & !anyNA(init) & is.vector(init) & length(init) == s ^ k)) {
    stop("'init' is not a numeric vector of length s ^ k")
  }
  
  if (!(all(init >= 0) & all(init <= 1))) {
    stop("Probabilities in 'init' must be between [0, 1]")
  }
  
  if (!((sum(init) >= 1 - sqrt(.Machine$double.eps)) | (sum(init) <= 1 + sqrt(.Machine$double.eps)))) {
    stop("The sum of 'init' is not equal to one")
  }
  
  #############################
  # Checking parameter k
  #############################
  
  if (!((k > 0) & ((k %% 1) == 0))) {
    stop("'k' must be a strictly positive integer")
  }
  
  #############################
  # Checking parameter ptrans
  #############################
  
  if (!(is.numeric(ptrans) & !anyNA(ptrans) & is.matrix(ptrans))) {
    stop("'ptrans' is not a matrix with numeric values")
  }
  
  if (!((dim(ptrans)[1] == s ^ k) & (dim(ptrans)[2] == s))) {
    stop("The dimension of the matrix 'ptrans' must be equal to (s ^ k, s)")
  }
  
  if (!(all(ptrans >= 0) & all(ptrans <= 1))) {
    stop("Probabilities in 'ptrans' must be between [0, 1]")
  }
  
  if (!all((apply(ptrans, 1, sum) >= 1 - sqrt(.Machine$double.eps)) | (apply(ptrans, 1, sum) <= 1 + sqrt(.Machine$double.eps)))) {
    stop("'ptrans' is not a stochastic matrix (column sums accross rows must be equal to one for each row)")
  }
  
  
  # Add names to the attributes init, ptrans, for readability
  colnames(ptrans) <- words(length = 1, alphabet = states)
  row.names(ptrans) <- words(length = k, alphabet = states)
  names(init) <- row.names(ptrans)
  
  ans <- list(states = states, s = s, k = k, init = init, ptrans = ptrans)
  
  class(ans) <- "mm"
  
  return(ans)
}


#' Function to check if an object is of class `mm`
#' 
#' @description `is.mm` returns `TRUE` if `x` is an object of 
#'   class `mm`.
#' 
#' @param x An arbitrary R object.
#' @return `is.mm` returns `TRUE` or `FALSE` depending on whether `x` is an 
#'   object of class `mm` or not.
#' 
#' @export
#' 
is.mm <- function(x) {
  inherits(x, "mm")
}


# Method to get the number of parameters
# (useful for the computation of criteria such as AIC and BIC)
.getKpar.mm <- function(x) {
  
  s <- x$s
  
  kpar <- (s - 1) * s ^ x$k
  
  return(kpar)
}


#' Log-likelihood Function
#' 
#' @description Computation of the log-likelihood for a Markov model
#' 
#' @param x An object of class [mm].
#' @param processes An object of class `processesMarkov`.
#' 
#' @noRd
#' 
.loglik.mm <- function(x, processes) {
  
  #############################
  # Let's compute the log-likelihood
  #############################
  
  Nstarti <- processes$Nstarti
  maskNstarti <- processes$Nstarti != 0 & x$init != 0
  
  Nij <- processes$Nij
  maskNij <- processes$Nij != 0 & x$ptrans != 0
  
  loglik <- sum(Nstarti[maskNstarti] * log(x$init[maskNstarti])) + sum(Nij[maskNij] * log(x$ptrans[maskNij]))
  
  return(loglik)
  
}


#' Akaike Information Criterion (AIC)
#' 
#' @description Computation of the Akaike Information Criterion.
#' 
#' @param x An object of class [mm].
#' @param sequences A list of vectors representing the sequences for which the 
#'   AIC will be computed based on `x`.
#' @return Value of the AIC.
#' 
#' @noRd
#' 
#' @export
#' 
aic.mm <- function(x, sequences) {
  
  loglik <- loglik(x, sequences)
  
  kpar <- .getKpar(x)
  
  aic <- -2 * loglik + 2 * kpar
  
  return(aic)
  
}


#' Bayesian Information Criterion (BIC)
#' 
#' @description Computation of the Bayesian Information Criterion.
#' 
#' @param x An object of class [mm].
#' @param sequences A list of vectors representing the sequences for which the 
#'   BIC will be computed based on `x`.
#' @return Value of the BIC.
#' 
#' @noRd
#' 
#' @export
#' 
bic.mm <- function(x, sequences) {
  
  loglik <- loglik(x, sequences)
  
  kpar <- .getKpar(x)
  
  n <- sum(sapply(sequences, length))
  
  bic <- -2 * loglik + log(n) * kpar
  
  return(bic)
}


#' Loglikelihood
#' 
#' @description Computation of the log-likelihood for a Markov model
#' 
#' @param x An object of class [mm].
#' @param sequences A list of vectors representing the sequences for which the 
#'   log-likelihood will be computed based on `x`.
#' @return Value of the log-likelihood.
#' 
#' @noRd
#' 
#' @export
#' 
loglik.mm <- function(x, sequences) {
  
  #############################
  # Checking parameters sequences and states
  #############################
  
  if (!(is.list(sequences) & all(sapply(sequences, class) %in% c("character", "numeric")))) {
    stop("The parameter 'sequences' should be a list of vectors")
  }
  
  if (!all(unique(unlist(sequences)) %in% x$states)) {
    stop("Some states in the list of observed sequences 'sequences' 
         are not in the state space given by the model 'x'")
  }
  
  processes <- processesMarkov(sequences = sequences, states = x$states, k = x$k, verbose = FALSE)
  loglik <- .loglik.mm(x = x, processes = processes)
  
  return(loglik)
  
}


#' Simulates k-th order Markov chains
#' 
#' @description Simulates k-th order Markov chains.
#' 
#' @details If `nsim` is a single integer then a chain of that length is 
#'   produced. If `nsim` is a vector of integers, then `length(nsim)` 
#'   sequences are generated with respective lengths.
#' 
#' @param object An object of class [mm].
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
#' @examples 
#' states <- c("a", "c", "g", "t")
#' s <- length(states)
#' k <- 2
#' init <- rep.int(1 / s ^ k, s ^ k)
#' p <- matrix(0.25, nrow = s ^ k, ncol = s)
#' 
#' # Specify a Markov model of order 1
#' markov <- mm(states = states, init = init, ptrans = p, k = k)
#' 
#' seqs <- simulate(object = markov, nsim = c(1000, 10000, 2000), seed = 150)
#' 
simulate.mm <- function(object, nsim = 1, seed = NULL, ...) {
  
  #############################
  # Checking parameter nsim
  #############################
  
  if (!all(is.numeric(nsim), is.vector(nsim), !anyNA(nsim), nsim > 0, (nsim %% 1) == 0)) {
    stop("'nsim' must be a strictly positive integer or a vector of striclty positive integers")
  }
  
  #############################
  # Checking parameter seed
  #############################
  
  if (is.null(seed)) {
    seed <- round(as.numeric(Sys.time()))
  }
  
  if (!all(is.numeric(seed), seed >= 0, (seed %% 1) == 0)) {
    stop("'seed' must be a positive integer")
  }
  
  
  s <- length(object$states)
  out <- list()
  nbseq <- length(nsim)
  
  for (n in 1:nbseq) {
    
    if (nsim[n] <= object$k) {
      
      y <- s2c(sample(x = names(object$init), size = 1, prob = object$init))
    
    } else {
    
      y <- rep.int(NA, nsim[n])
      
      # Initial state(s)
      y[1:object$k] <- s2c(sample(x = names(object$init), size = 1, prob = object$init))
      
      for (i in 1:(nsim[n] - object$k)) {
        ind <- which(object$states == y[i + object$k - 1])
        if (object$k > 1) {
          for (j in (object$k - 2):0) {
            ind <- ind + s ^ (j + 1) * (which(object$states == y[i + j]) - 1)
          }
        }
        y[i + object$k] <- sample(object$states, 1, prob = object$ptrans[ind, ])
      }
      
    }
    
    out[[n]] <- y
    
  }
  
  return(out)
  
}
