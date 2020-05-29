#' Estimation of a k-th order Markov chain
#'
#' @description Estimation of the transition matrix and initial law of a k-th 
#'   order Markov chain starting from one or several sequences.
#'
#' @details Let \eqn{X_1, X_2, ..., X_n} be a trajectory of length \eqn{n} of 
#'   the Markov chain \eqn{X = (X_m)_{m \in N}} of order \eqn{k = 1} with 
#'   transition matrix \eqn{p_{trans}(i,j) = P(X_{m+1} = j | X_m = i)}. The 
#'   estimation of the transition matrix is \eqn{\widehat{p_{trans}}(i,j) = \frac{N_{ij}}{N_{i.}}}, 
#'   where \eqn{N_{ij}} is the number of transitions from state \eqn{i} to state 
#'   \eqn{j} and \eqn{N_{i.}} is the number of transition from state \eqn{i} 
#'   to any state. For \eqn{k > 1} we have similar expressions.
#'
#'  The initial distribution of a k-th order Markov chain is defined as 
#'  \eqn{\mu_i = P(X_1 = i)}. An estimation of the initial law for a first order 
#'  Markov chain is assumed to be the estimation of the stationary distribution. 
#'  If the order of the Markov is greater than 1, then an estimation of the 
#'  initial law is \eqn{\widehat{\mu_i} = \frac{N_i}{N}}, where \eqn{N_i} is the number occurences 
#'  of state \eqn{i} in the sequences and \eqn{N} is the sum of the sequence 
#'  lengths.
#'
#' @param sequences A list of vectors representing the sequences.
#' @param states Vector of state space (of length s).
#' @param k Order of the Markov chain.
#' @return An object of class [markovmodel][markovmodel].
#' 
#' 
#' @seealso [markovmodel], [simulate], [smmnonparametric], [smmparametric], [fitsemimarkovmodel]
#' @export
#'
#' @examples 
#' states <- c("a", "c", "g", "t")
#' s <- length(states)
#' vect.init <- c(1 / 4, 1 / 4, 1 / 4, 1 / 4)
#' k <- 2
#' p <- matrix(0.25, nrow = s ^ k, ncol = s)
#' 
#' # Specify the Markov model
#' markov1 <- markovmodel(states = states, init = vect.init, ptrans = p, k = k)
#' 
#' seq1 <- simulate(object = markov1, nsim = c(1000, 10000, 2000), seed = 150)
#' 
#' est <- fitmarkovmodel(sequences = seq1, states = states, k = 2)
#' 
fitmarkovmodel <- function(sequences, states, k = 1) {
  
  #############################
  # Checking parameters sequences and states
  #############################
  
  if (!(is.list(sequences))) {
    stop("The parameter sequences should be a list")
  }
  
  if (!all(unique(unlist(sequences)) %in% states)) {
    stop("Some states in the list of observed sequences sequences are not in the state space states")
  }
  
  #############################
  # Checking parameter k
  #############################
  
  if (!((k > 0) && ((k %% 1) == 0))) {
    stop("k must be a strictly positive integer")
  }
  
  
  s <- length(states) # State space size
  nbseq <- length(sequences) # Number of sequences
  sequence <- c()
  
  # Count the number of transitions from i to j in k steps
  Nijl <- array(0, c(s ^ k, s, nbseq))
  
  for (i in 1:nbseq) {
    
    sequence <- sequences[[i]]
    Nijl[, , i] <- matrix(count(seq = sequence, wordsize = k + 1, alphabet = states), byrow = TRUE, ncol = s)
    
  }
  
  Nij <- apply(Nijl, c(1, 2), sum)
  Ni <- apply(Nijl, 1, sum)
  
  indexdiag <- seq(1, s * s, by = s + 1)
  Nij.temp <- as.vector(Nij)[-indexdiag]
  
  statesi <- row(Nij)[-indexdiag][which(Nij.temp == 0)]
  statesj <- col(Nij)[-indexdiag][which(Nij.temp == 0)]
  
  # Verify the existence of all transitions
  if ((length(statesi) != 0) | (length(statesj) != 0)) {
    warning(
      "Some transitions from state i to state j are not observed.
            The following are ",
      paste0(sapply(1:length(statesi), function(x)
        paste0("(i=", statesi[x], " to j=", statesj[x], ")")), collapse = ", "),
      "."
    )
  }
  
  # Verify the existence of all states
  if (length(which(Ni == 0)) > 0) {
    warning("Missing observed states")
  }
  
  # Compute the transition matrix
  ptrans <- Nij / tcrossprod(Ni, rep.int(1, s))
  ptrans[which(is.na(ptrans))] <- 0
  
  if (k == 1) {
    init <- .stationaryDistribution(ptrans)
  } else {
    Nstart <- as.vector(count(seq = unlist(sequences), wordsize = 1, alphabet = states)) ## initial law for k > 1 ???
    init <- Nstart / sum(Nstart)
  }
  
  estimate <- list(states = states, s = s, init = init, ptrans = ptrans, k = k)
  class(estimate) <- "markovmodel"
  
  return(estimate)
}
