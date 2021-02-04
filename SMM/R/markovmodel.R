#' Markov model specification
#'
#' @description Creates a model specification of a Markov model.
#'
#' @param states Vector of state space of length s.
#' @param k Order of the Markov chain.
#' @param init Vector of initial distribution of length s ^ k.
#' @param ptrans Matrix of transition probabilities of dimension \eqn{(s, s)}.
#' @return An object of class [markovmodel][markovmodel].
#' 
#' @seealso [simulate.markovmodel], [fitmarkovmodel]
#' 
#' @export
#'
markovmodel <- function(states, init, ptrans, k = 1) {

  #############################
  # Checking parameter states
  #############################

  s <- length(states)
  
  if (!(is.vector(states) && (length(unique(states)) == s))) {
    stop("The state space states is not a vector of unique elements")
  }

  #############################
  # Checking parameter init
  #############################

  if (!(is.vector(init) && (length(init) == s ^ k))) {
    stop("init is not a vector of length s ^ k.")
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
  
  if (!((dim(ptrans)[1] == s ^ k) && (dim(ptrans)[2] == s))) {
    stop("The dimension of the matrix ptrans must be equal to (s ^ k, s)")
  }
  
  if (!all(apply(ptrans, 1, sum) == 1)) {
    stop("ptrans is not a stochastic matrix (column sums accross rows must be equal to one for each row)")
  }
  
  
  # Add names to the attributes init, ptrans, for readability
  colnames(ptrans) <- words(length = 1, alphabet = states)
  row.names(ptrans) <- words(length = k, alphabet = states)
  names(init) <- row.names(ptrans)

  ans <- list(states = states, s = s, k = k, init = init, ptrans = ptrans)
  
  class(ans) <- "markovmodel"
  
  return(ans)
}

# Function to check if an object is of class markovmodel
is.markovmodel <- function(x) {
  inherits(x, "markovmodel")
}

#' Loglikelihood
#'
#' @description Computation of the log-likelihood for a Markov model
#'
#' @param x An object of class [markovmodel].
#' @param sequences A list of vectors representing the sequences for which the 
#'   log-likelihood will be computed based on `x`.
#' @return Value of the log-likelihood.
#' 
#' @export
#'
loglik.markovmodel <- function(x, sequences) {
  
  #############################
  # Checking parameters sequences and states
  #############################
  
  if (!(is.list(sequences) && all(sapply(sequences, class) %in% c("character", "numeric")))) {
    stop("The parameter sequences should be a list of vectors")
  }
  
  if (!all(unique(unlist(sequences)) %in% x$states)) {
    stop("Some states in the list of observed sequences sequences are not in the state space given by the model x")
  }
  
  s <- length(x$states)
  nbseq <- length(sequences) # Number of sequences
  
  # Count the number of transitions from i to j in k steps
  Nijl <- array(0, c(s ^ x$k, s, nbseq))
  contrInit <- rep.int(x = 0, times = nbseq)
  
  for (i in 1:nbseq) {
    
    # Transitions
    Nijl[, , i] <- matrix(count(seq = sequences[[i]], wordsize = x$k + 1, alphabet = x$states), byrow = TRUE, ncol = s)
    
    # Initial state(s)
    contrInit[i] <- log(x$init[which(words(length = x$k, alphabet = x$states) == c2s(sequences[[i]][1:x$k]))])
  }
  
  Nij <- apply(Nijl, c(1, 2), sum)
  maskNij <- Nij != 0 & x$ptrans != 0
  
  loglik <- sum(contrInit) + sum(Nij[maskNij] * log(x$ptrans[maskNij]))
  
  return(loglik)
}

# Method to get the number of parameters
# (useful for the computation of criteria such as AIC and BIC)
.getKpar.markovmodel <- function(x) {
  
  s <- x$s
  
  kpar <- (s - 1) * s ^ x$k
  
  return(kpar)
}

#' Akaike Information Criterion (AIC)
#'
#' @description Computation of the Akaike Information Criterion.
#'
#' @param x An object of class [markovmodel].
#' @param sequences A list of vectors representing the sequences for which the 
#'   AIC will be computed based on `x`.
#' @return Value of the AIC.
#' 
#' @export
#'
aic.markovmodel <- function(x, sequences) {
  
  loglik <- loglik(x, sequences)
  
  kpar <- .getKpar(x)
  
  aic <- -2 * loglik + 2 * kpar
  
  return(aic)
  
}

#' Bayesian Information Criterion (BIC)
#'
#' @description Computation of the Bayesian Information Criterion.
#'
#' @param x An object of class [markovmodel].
#' @param sequences A list of vectors representing the sequences for which the 
#'   BIC will be computed based on `x`.
#' @return Value of the BIC.
#' 
#' @export
#'
bic.markovmodel <- function(x, sequences) {
  
  loglik <- loglik(x, sequences)
  
  kpar <- .getKpar(x)
  
  n <- sum(unlist(lapply(sequences, length)))
  
  bic <- -2 * loglik + log(n) * kpar
  
  return(bic)
}
