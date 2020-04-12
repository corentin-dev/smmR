is.smm <- function(x) {
  inherits(x, "smm")
}

# Method to get the sojourn time distribution f
.get.f <- function(x, ...) {
  UseMethod(".get.f", x)
}

# Method to get the kernel q
.get.q <- function(x, ...) {
  UseMethod(".get.q", x)
}

#' Loglikelihood
#'
#' @description Computation of the loglikelihood for a semi-Markov model
#'
#' @param x An object inheriting from the class smm
#'   ([smmnonparametric][smmnonparametric] or [smmparametric][smmparametric]).
#' @param seq A list of vectors representing the sequences for which the 
#'   log-likelihood must be computed.
#' @param E Vector of state space (of length S).
#' @return A vector giving the value of the loglikelihood for each sequence.
#' 
#' 
#' @export
#'
loglik.smm <- function(x, seq, E) {
  
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
  # Checking smm parameter
  #############################
  
  if ((x$S != S)) {
    stop("The size of the matrix ptrans must be equal to SxS with S = length(E)")
  }
  
  if (!all.equal(E, x$E)) {
    stop("The state space of the estimated SMM smm is different from the given state E")
  }
  
  
  seq <- sequences(listSeq = seq, E = E)
  Kmax <- seq$Kmax
  
  if (!(is.null(x$Kmax))) {
    if (!(Kmax == x$Kmax)) {
      stop("Kmax of the given sequences is different from the Kmax of the estimated SMM model")  
    }
  }
  
  type.sojourn <- x$type.sojourn
  init <- x$init
  pij <- x$ptrans
  
  fijk <- .get.f(x, Kmax = Kmax)
  
  loglik <- rep.int(x = NA, times = seq$nbSeq)
  
  for (l in 1:seq$nbSeq) {
    Nij <- seq$counting$Nijl[, , l]
    Nstart <- seq$counting$Nstart[, l]
    Nijk <- seq$counting$Nijkl[, , , l]
    
    mask <- Nijk != 0 && fijk != 0
    
    loglik[l] <- sum(Nstart * log(init)) +
      sum(Nij[row(Nij) != col(Nij)] * log(pij[row(pij) != col(pij)])) +
      sum(Nijk[mask] * log(fijk[mask]))
    
  }
  
  return(loglik)
}

# Method to get the number of parameters
.getKpar <- function(x) {
  UseMethod(".getKpar", x)
}

#' Akaike Information Criterion (AIC)
#'
#' @description Computation of the Akaike Information Criterion.
#'
#' @param x An object inheriting from the class smm
#'   ([smmnonparametric][smmnonparametric] or [smmparametric][smmparametric]).
#' @param seq A list of vectors representing the sequences for which the 
#'   AIC criterion must be computed.
#' @param E Vector of state space (of length S).
#' @return A vector giving the value of the AIC for each sequence.
#' 
#' 
#' @export
#'
aic.smm <- function(x, seq, E) {
  
  loglik <- loglik(x, seq, E)
  seq <- sequences(listSeq = seq, E = E)
  
  S <- x$S
  Kmax <- seq$Kmax
  
  Kpar <- .getKpar(x)
  
  aic <- rep.int(x = NA, times = seq$nbSeq)
  for (l in 1:seq$nbSeq) {
    aic[l] <- -2 * loglik[l] + 2 * Kpar
  }
  
  return(aic)
  
}

#' Bayesian Information Criterion (BIC)
#'
#' @description Computation of the Bayesian Information Criterion.
#'
#' @param x An object inheriting from the class smm
#'   ([smmnonparametric][smmnonparametric] or [smmparametric][smmparametric]).
#' @param seq A list of vectors representing the sequences for which the 
#'   BIC criterion must be computed.
#' @param E Vector of state space (of length S).
#' @return A vector giving the value of the BIC for each sequence.
#' 
#' 
#' @export
#'
bic.smm <- function(x, seq, E) {
  
  loglik <- loglik(x, seq, E)
  seq <- sequences(listSeq = seq, E = E)
    
  S <- x$S
  Kmax <- seq$Kmax
    
  Kpar <- .getKpar(x)
  
  bic <- rep.int(x = NA, times = seq$nbSeq)
  for (l in 1:seq$nbSeq) {
    n <- length(seq$J[[l]])
    bic[l] <- -2 * loglik[l] + log(n) * Kpar
  }
  
  return(bic)
    
}