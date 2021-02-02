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
#'  \eqn{\mu_i = P(X_1 = i)}. Three methods are proposed for the estimation
#'  of the latter :
#'  \describe{
#'    \item{Estimation based on the Maximum Likelihood Estimator:}{
#'      The Maximum Likelihood Estimator for the initial distribution. The 
#'      formula is: \eqn{\widehat{\mu_i} = \frac{Nstart_i}{L}}, where 
#'      \eqn{Nstart_i} is the number of occurences of the word \eqn{i} (of 
#'      length \eqn{k}) at the beginning of each sequence and \eqn{L} is the 
#'      number of sequences. This estimator is reliable when the number of 
#'      sequences \eqn{L} is high.}
#'    \item{Estimation based on the frequency:}{The initial distribution is 
#'      estimated by taking the frequences of the words of length `k` for all 
#'      sequences. The formula is \eqn{\widehat{\mu_i} = \frac{N_i}{N}}, where 
#'      \eqn{N_i} is the number of occurences of the word \eqn{i} (of length \eqn{k}) 
#'      in the sequences and \eqn{N} is the sum of the lengths of the sequences.}
#'    \item{Estimation based on the product of the frequences of each state:}{
#'      The initial distribution is estimated by using the product of the 
#'      frequences of each state (for all the sequences) in the word of length
#'      \eqn{k}.}
#'  }
#'
#' @param sequences A list of vectors representing the sequences.
#' @param states Vector of state space (of length s).
#' @param k Order of the Markov chain.
#' @param init.estim Optional. Method used to estimate the initial distribution.
#'   If `init.estim = "mle"`, then the classical Maximum Likelihood Estimator 
#'   is used, if `init.estim = "freq"`, then, the initial distribution `init` 
#'   is estimated by taking the frequences of the words of length `k` for all 
#'   sequences. If `init.estim = "prod"`, then, `init` is estimated by using 
#'   the product of the frequences of each letter (for all the sequences) in 
#'   the word of length `k`. See Details for the formulas.
#' @return An object of class [markovmodel][markovmodel].
#' 
#' 
#' @seealso [markovmodel], [simulate], [smmnonparametric], [smmparametric], [fitsemimarkovmodel]
#' @export
#'
#' @examples 
#' states <- c("a", "c", "g", "t")
#' s <- length(states)
#' k <- 2
#' vect.init <- rep.int(1 / s ^ k, s ^ k)
#' p <- matrix(0.25, nrow = s ^ k, ncol = s)
#' 
#' # Specify the Markov model
#' markov1 <- markovmodel(states = states, init = vect.init, ptrans = p, k = k)
#' 
#' seq1 <- simulate(object = markov1, nsim = c(1000, 10000, 2000), seed = 150)
#' 
#' est <- fitmarkovmodel(sequences = seq1, states = states, k = 2)
#' 
fitmarkovmodel <- function(sequences, states, k = 1, init.estim = c("mle", "freq", "prod")) {
  
  #############################
  # Checking parameters sequences and states
  #############################
  
  if (!(is.list(sequences) && all(sapply(sequences, class) %in% c("character", "numeric")))) {
    stop("The parameter sequences should be a list of vectors")
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
  
  #############################
  # Checking parameter init.estim
  #############################
  
  init.estim <- match.arg(init.estim)

  
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
  
  if (init.estim == "mle") {
    Nstart <- as.vector(count(seq = unlist(lapply(sequences, function(x) x[1:k])), wordsize = k, by = k, alphabet = states))
    init <- Nstart / sum(Nstart)
  } else if (init.estim == "freq") {
    Nstart <- as.vector(count(seq = unlist(sequences), wordsize = k, alphabet = states))
    init <- Nstart / sum(Nstart)
  } else {# init.estim == "prod"
    Nstart <- as.vector(count(seq = unlist(sequences), wordsize = 1, alphabet = states))
    prob <- Nstart / sum(Nstart)
    
    init <- as.vector(.productProb(length = k, prob = prob))
  }
  
  estimate <- markovmodel(states = states, init = init, ptrans = ptrans, k = k)
  
  return(estimate)
}
