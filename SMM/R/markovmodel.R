#' Markov model specification
#'
#' @description Creates a model specification of a Markov model.
#'
#' @param E Vector of state space of length S.
#' @param init Vector of initial distribution of length S.
#' @param ptrans Matrix of transition probabilities of the embedded Markov chain 
#'   \eqn{J=(J_m)_{m}} of size SxS.
#' @param k Order of the Markov chain.
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


  ans <- list(E = E, init = init, ptrans = ptrans, k = k)
  
  class(ans) <- "markovmodel"
  
  return(ans)
}

is.markovmodel <- function(x) {
  inherits(x, "markovmodel")
}

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
  
  for (i in 1:nbseq) {
    
    Nijl[, , i] <- matrix(count(seq = seq[[i]], wordsize = x$k + 1, alphabet = x$E), byrow = TRUE, ncol = x$S)
    
  }
  
  loglik <- rep.int(x = NA, times = nbseq)
  for (j in 1:nbseq) {
    s <- 0
    for (i in 1:x$k) {# Warning to initial law
      s <- s + log(x$init[which(x$E == seq[[j]][i])])
    }
    loglik[j] <- s + sum(as.numeric(Nijl[, , j])[which(x$ptrans != 0)] * log(x$ptrans[which(x$ptrans != 0)]))
  }
  
  return(loglik)
}

# Method to get the number of parameters
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
#' @return A vector giving the value of the AIC for each sequence.
#' 
#' 
#' @export
#'
aic.markovmodel <- function(x, seq, E) {
  
  loglik <- loglik(x, seq, E)
  nbseq <- length(seq)
  
  Kpar <- .getKpar(x)
  
  aic <- rep.int(x = NA, times = nbseq)
  for (l in 1:nbseq) {
    aic[l] <- -2 * loglik[l] + 2 * Kpar
  }
  
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
#' @return A vector giving the value of the BIC for each sequence.
#' 
#' 
#' @export
#'
#' @examples
#' E <- c("a", "c", "g", "t")
#' S <- length(E)
#' vect.init <- c(1 / 4, 1 / 4, 1 / 4, 1 / 4)
#' k <- 2
#' p <- matrix(0.25, nrow = S ^ k, ncol = S)
#' 
#' # Specify the Markov model
#' markov1 <- markovmodel(E = E, init = vect.init, ptrans = p, k = k)
#' markov1
#' 
bic.markovmodel <- function(x, seq, E) {
  
  out <- loglik(x, seq, E)
  nbseq <- length(seq)
  loglik <- out$loglik
  
  Kpar <- .getKpar(x)
  
  bic <- rep.int(x = NA, times = nbseq)
  for (l in 1:nbseq) {
    n <- length(seq[[l]])
    bic[l] <- -2 * loglik[l] + log(n) * Kpar
  }
  
  return(bic)
}
