#' Markov model specification
#'
#' @description Creates a model specification of a Markov model.
#'
#' @param E Vector of state space of length S.
#' @param k Order of the Markov chain.
#' @param init Vector of initial distribution of length S.
#' @param ptrans Matrix of transition probabilities of the embedded Markov chain 
#'   \eqn{J=(J_m)_{m}} of size SxS.
#' @return An object of class [markovmodel][markovmodel].
#' 
#' 
#' @seealso [simulate], [fitmarkovmodel], [smmnonparametric], [smmparametric], [fitsemimarkovmodel]
#' @export
#'
markovmodel <- function(E, init, ptrans, k = 1) {

  #############################
  # Checking parameter E
  #############################

  S <- length(E)
  
  if (!(is.vector(E) && (length(unique(E)) == S))) {
    stop("The state space E is not a vector of unique elements")
  }

  #############################
  # Checking parameter init
  #############################

  if (!(is.vector(init) && (length(init) == S))) {
    stop("init is not a vector of length S")
  }
  
  if (!(all(init >= 0) && all(init <= 1))) {
    stop("Probabilities in init must be between [0, 1]")
  }
  
  if (!(sum(init) == 1)) {
    stop("The sum of init is not equal to one")
  }

  #############################
  # Checking parameter k
  #############################
  
  if (!((k > 0) && ((k %% 1) == 0))) {
    stop("k must be a strictly positive integer")
  }
  
  #############################
  # Checking parameter ptrans
  #############################

  if (!is.matrix(ptrans)) {
    stop("ptrans is not a matrix")
  }
  
  if (!(all(ptrans >= 0) && all(ptrans <= 1))) {
    stop("Probabilities in ptrans must be between [0, 1]")
  }
  
  if (!((dim(ptrans)[1] == S ^ k) && (dim(ptrans)[2] == S))) {
    stop("The size of the matrix ptrans must be equal to SxS")
  }
  
  if (!all(apply(ptrans, 1, sum) == 1)) {
    stop("ptrans is not a stochastic matrix (column sums accross rows must be equal to one for each row)")
  }


  ans <- list(E = E, k = k, init = init, ptrans = ptrans)
  
  class(ans) <- "markovmodel"
  
  return(ans)
}

# Function to check if an object is of class markovmodel
is.markovmodel <- function(x) {
  inherits(x, "markovmodel")
}

# Method to compute the stationary distribution of the Markov model x
.stationary.distribution.markovmodel <- function(x) {
  
  m <- dim(x$ptrans)[1] # Number of states
  
  A <- t(x$ptrans) - diag(1, m, m)
  A[m, ] <- 1
  b <- c(rep(0, (m - 1)), 1)
  statlaw <- solve(A, b)
  
  return(statlaw)
}

#' Loglikelihood
#'
#' @description Computation of the loglikelihood for a semi-Markov model
#'
#' @param x An object of class [markovmodel][markovmodel].
#' @param seq A list of vectors representing the sequences for which the 
#'   log-likelihood must be computed.
#' @param E Vector of state space (of length S).
#' @return A vector giving the value of the loglikelihood for each sequence.
#' 
#' 
#' @export
#'
loglik.markovmodel <- function(x, seq, E) {
  
  #############################
  # Checking parameters seq and E
  #############################
  
  if (!is.list(seq)) {
    stop("The parameter seq should be a list")
  }
  
  if (!all(unique(unlist(seq)) %in% E)) {
    stop("Some states in the list of observed sequences seq are not in the state space E")
  }
  
  S <- length(E)
  
  #############################
  # Checking markovmodel parameter
  #############################
  
  if ((x$S != S)) {
    stop("The size of the matrix ptrans must be equal to SxS with S = length(E)")
  }
  
  if (!all.equal(E, x$E)) {
    stop("The state space of the estimated Markov model is different from the given state E")
  }
  
  
  nbseq <- length(seq) # Number of sequences
  vect.seq <- c()
  
  # Count the number of transitions from i to j in k steps
  Nijl <- array(0, c(S ^ x$k, S, nbseq))
  contrInit <- rep.int(x = 0, times = nbseq)
  
  for (i in 1:nbseq) {
    
    Nijl[, , i] <- matrix(count(seq = seq[[i]], wordsize = x$k + 1, alphabet = x$E), byrow = TRUE, ncol = x$S)
    
    for (j in 1:x$k) {# Warning to initial law
      if (x$init[which(x$E == seq[[i]][j])] != 0) {
        contrInit[i] <- contrInit[i] + log(x$init[which(x$E == seq[[i]][j])])  
      }
    }
    
  }
  
  Nij <- apply(Nijl, c(1, 2), sum)
  maskNij <- Nij != 0 & x$ptrans != 0
  
  loglik <- sum(contrInit) + sum(Nij[maskNij] * log(x$ptrans[maskNij]))
  
  return(loglik)
}

# Method to get the number of parameters
# (useful for the computation of criteria such as AIC and BIC)
.getKpar.markovmodel <- function(x) {
  
  Kpar <- (x$S - 1) * x$S ^ x$k
  
  return(Kpar)
}

#' Akaike Information Criterion (AIC)
#'
#' @description Computation of the Akaike Information Criterion.
#'
#' @param x An object of class [markovmodel][markovmodel].
#' @param seq A list of vectors representing the sequences for which the 
#'   AIC criterion must be computed.
#' @param E Vector of state space (of length S).
#' @return A numeric value giving the value of the AIC.
#' 
#' 
#' @export
#'
aic.markovmodel <- function(x, seq, E) {
  
  loglik <- loglik(x, seq, E)
  
  Kpar <- .getKpar(x)
  
  aic <- -2 * loglik + 2 * Kpar
  
  return(aic)
  
}

#' Bayesian Information Criterion (BIC)
#'
#' @description Computation of the Bayesian Information Criterion.
#'
#' @param x An object of class [markovmodel][markovmodel].
#' @param seq A list of vectors representing the sequences for which the 
#'   BIC criterion must be computed.
#' @param E Vector of state space (of length S).
#' @return A numeric value giving the value of the BIC.
#' 
#' 
#' @export
#'
bic.markovmodel <- function(x, seq, E) {
  
  loglik <- loglik(x, seq, E)
  
  Kpar <- .getKpar(x)
  
  n <- sum(unlist(lapply(seq, length)))
  bic <- -2 * loglik + log(n) * Kpar
  
  return(bic)
}
