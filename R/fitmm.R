#' Maximum Likelihood Estimation (MLE) of a k-th order Markov chain
#' 
#' @description Maximum Likelihood Estimation of the transition matrix and 
#'   initial distribution of a k-th order Markov chain starting from one or 
#'   several sequences.
#' 
#' @details Let \eqn{X_1, X_2, ..., X_n} be a trajectory of length \eqn{n} of 
#'   the Markov chain \eqn{X = (X_m)_{m \in N}} of order \eqn{k = 1} with 
#'   transition matrix \eqn{p_{trans}(i,j) = P(X_{m+1} = j | X_m = i)}. The 
#'   maximum likelihood estimation of the transition matrix is 
#'   \eqn{\widehat{p_{trans}}(i,j) = \frac{N_{ij}}{N_{i.}}}, where \eqn{N_{ij}}
#'   is the number of transitions from state \eqn{i} to state \eqn{j} and 
#'   \eqn{N_{i.}} is the number of transition from state \eqn{i} to any state. 
#'   For \eqn{k > 1} we have similar expressions.
#'   
#'  The initial distribution of a k-th order Markov chain is defined as 
#'  \eqn{\mu_i = P(X_1 = i)}. Five methods are proposed for the estimation
#'  of the latter :
#'  \describe{
#'    \item{Maximum Likelihood Estimator: }{The Maximum Likelihood Estimator 
#'      for the initial distribution. The formula is: 
#'      \eqn{\widehat{\mu_i} = \frac{Nstart_i}{L}}, where \eqn{Nstart_i} is 
#'      the number of occurences of the word \eqn{i} (of length \eqn{k}) at 
#'      the beginning of each sequence and \eqn{L} is the number of sequences. 
#'      This estimator is reliable when the number of sequences \eqn{L} is high.}
#'    \item{Stationary distribution: }{The stationary distribution is used as 
#'      a surrogate of the initial distribution. If the order of the Markov 
#'      chain is more than one, the transition matrix is converted into a 
#'      square block matrix in order to estimate the stationary distribution.
#'      This method may take time if the order of the Markov chain is high 
#'      (more than three (3)).}
#'    \item{Frequencies of each word: }{The initial distribution is estimated 
#'      by taking the frequencies of the words of length `k` for all sequences.
#'      The formula is \eqn{\widehat{\mu_i} = \frac{N_i}{N}}, where \eqn{N_i} 
#'      is the number of occurences of the word \eqn{i} (of length \eqn{k}) in 
#'      the sequences and \eqn{N} is the sum of the lengths of the sequences.}
#'    \item{Product of the frequencies of each state: }{The initial distribution 
#'      is estimated by using the product of the frequencies of each state 
#'      (for all the sequences) in the word of length \eqn{k}.}
#'    \item{Uniform distribution: }{The initial probability of each state is 
#'       equal to \eqn{1 / s}, with \eqn{s}, the number of states.}
#'  }
#'  
#' @param sequences A list of vectors representing the sequences.
#' @param states Vector of state space (of length s).
#' @param k Order of the Markov chain.
#' @param init.estim Optional. `init.estim` gives the method used to estimate 
#'   the initial distribution. The following methods are proposed:
#'   \itemize{
#'     \item `init.estim = "mle"`: the classical Maximum Likelihood Estimator 
#'       is used to estimate the initial distribution `init`;
#'     \item `init.estim = "stationary"`: the initial distribution is replaced by 
#'       the stationary distribution of the Markov chain (if the order of the 
#'       Markov chain is more than one, the transition matrix is converted 
#'       into a square block matrix in order to estimate the stationary 
#'       distribution);
#'     \item `init.estim = "freq"`: the initial distribution is estimated by 
#'       taking the frequencies of the words of length `k` for all sequences;
#'     \item `init.estim = "prod"`: `init` is estimated by using the product 
#'       of the frequencies of each letter (for all the sequences) in the word 
#'       of length `k`;
#'     \item `init.estim = "unif"`: the initial probability of each state is 
#'       equal to \eqn{1 / s}, with \eqn{s} the number of states.
#'   }
#'   
#' @return An object of class S3 `mmfit` (inheriting from the S3 class [mm]).
#'   The S3 class `mmfit` contains:
#'   \itemize{
#'     \item All the attributes of the S3 class [mm];
#'     \item An attribute `M` which is an integer giving the total length of 
#'       the set of sequences `sequences` (sum of all the lengths of the list 
#'       `sequences`);
#'     \item An attribute `loglik` which is a numeric value giving the value 
#'       of the log-likelihood of the specified Markov model based on the 
#'       `sequences`;
#'     \item An attribute `sequences` which is equal to the parameter 
#'       `sequences` of the function `fitmm` (i.e. a list of sequences used to 
#'       estimate the Markov model).
#'   }
#'   
#' 
#' @seealso [mm], [simulate.mm]
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
#' # Specify a Markov model of order 2
#' markov <- mm(states = states, init = init, ptrans = p, k = k)
#' 
#' seqs <- simulate(object = markov, nsim = c(1000, 10000, 2000), seed = 150)
#' 
#' est <- fitmm(sequences = seqs, states = states, k = 2)
#' 
fitmm <- function(sequences, states, k = 1, init.estim = "mle") {
  
  #############################
  # Checking parameters sequences and states
  #############################
  
  if (!(is.list(sequences) & all(sapply(sequences, class) %in% c("character", "numeric")))) {
    stop("The parameter 'sequences' should be a list of vectors")
  }
  
  if (!all(unique(unlist(sequences)) %in% states)) {
    stop("Some states in the list of observed sequences 'sequences' are not in the state space 'states'")
  }
  
  #############################
  # Checking parameter k
  #############################
  
  if (!((k > 0) & ((k %% 1) == 0))) {
    stop("'k' must be a strictly positive integer")
  }
  
  #############################
  # Checking parameter init.estim
  #############################
  
  # init.estim <- match.arg(init.estim)
  
  processes <- processesMarkov(sequences = sequences, states = states, k = k)
  s <- processes$s
  Nij <- processes$Nij
  Ni <- processes$Ni
  Nstarti <- processes$Nstarti
  
  # Compute the transition matrix
  ptrans <- Nij / tcrossprod(Ni, rep.int(1, s))
  ptrans[which(is.na(ptrans))] <- 0
  ptrans <- .normalizePtrans(ptrans)
  
  # Initial distribution
  if (is.vector(init.estim) & length(init.estim) == 1) {
    if (init.estim == "mle") {
      init <- Nstarti / sum(Nstarti)
    } else if (init.estim == "stationary") {
      if (k == 1) {
        init <- .stationaryDistribution(ptrans = ptrans)
      } else {
        init <- .stationaryDistribution(ptrans = .blockMatrix(ptrans = ptrans))
      }
    } else if (init.estim == "freq") {
      Nstart <- as.vector(count(seq = unlist(sequences), wordsize = k, alphabet = states))
      init <- Nstart / sum(Nstart)
    } else if (init.estim == "prod") {
      Nstart <- as.vector(count(seq = unlist(sequences), wordsize = 1, alphabet = states))
      prob <- Nstart / sum(Nstart)
      init <- as.vector(.productProb(length = k, prob = prob))
    } else if (init.estim == "unif") {
      init <- rep.int(x = 1 / (s ^ k), times = s ^ k)
    } else {
      stop("'init.estim' must be equal to \"mle\", \"stationary\", \"freq\", \"prod\" or \"unif\".
           'init.estim' can also be a vector of length s ^ k for custom initial distribution")
    }
  } else {
    if (!(is.numeric(init.estim) & !anyNA(init.estim) & is.vector(init.estim) & length(init.estim) == s ^ k)) {
      stop("'init.estim' is not a numeric vector of length s ^ k")
    }
    
    if (!(all(init.estim >= 0) & all(init.estim <= 1))) {
      stop("Probabilities in 'init.estim' must be between [0, 1]")
    }
    
    if (!((sum(init.estim) >= 1 - sqrt(.Machine$double.eps)) | (sum(init.estim) <= 1 + sqrt(.Machine$double.eps)))) {
      stop("The sum of 'init.estim' is not equal to one")
    }
    
    init <- init.estim
  }
  
  init <- as.vector(init / sum(init))
  
  mm <- mm(states = states, init = init, ptrans = ptrans, k = k)
  
  if (any(mm$init == 0)) {
    message("The probabilities of the initial state(s) \"", 
            paste0(names(which(mm$init == 0)), collapse = "\", \""),
            "\" are 0.")
  }
  
  loglik <- .loglik(x = mm, processes = processes)
  estimate <- mmfit(mm = mm, M = processes$M, loglik = loglik, sequences = sequences)
  
  return(estimate)
}
