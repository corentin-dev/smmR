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

# Method to get the log-likelihood
.loglik <- function(x, ...) {
  UseMethod(".loglik", x)
}

.loglik.smm <- function(x, listSeq, E) {
  
  #############################
  # Checking parameters listSeq and E
  #############################
  
  if (!is.list(listSeq)) {
    stop("The parameter listSeq should be a list")
  }
  
  if (!all(unique(unlist(listSeq)) %in% E)) {
    stop("Some states in the list of observed sequences listSeq are not in the state space E")
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
  
  
  seq <- sequences(listSeq = listSeq, E = E)
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
  
  return(list(seq = seq, loglik = loglik))
}

# Method to get the number of parameters
.getKpar <- function(x) {
  UseMethod(".getKpar", x)
}

# Method to get the AIC
AICSMM <- function(x, listSeq, E) {
  UseMethod("AICSMM", x)
}

AICSMM.smm <- function(x, listSeq, E) {
  
  out <- .loglik(x, listSeq, E)
  seq <- out$seq
  loglik <- out$loglik
  
  S <- x$S
  Kmax <- seq$Kmax
  
  Kpar <- .getKpar(x)
  
  aic <- rep.int(x = NA, times = seq$nbSeq)
  for (l in 1:seq$nbSeq) {
    aic[l] <- -2 * loglik[l] + 2 * Kpar
  }
  
  return(aic)
  
}

# Method to get the AIC
BICSMM <- function(x, listSeq, E) {
  UseMethod("BICSMM", x)
}

BICSMM.smm <- function(x, listSeq, E) {
  
  out <- .loglik(x, listSeq, E)
  seq <- out$seq
  loglik <- out$loglik
    
  S <- x$S
  Kmax <- seq$Kmax
    
  Kpar <- .getKpar(x)
  
  bic <- rep.int(x = NA, times = seq$nbSeq)
  for (l in 1:seq$nbSeq) {
    n <- length(listSeq[[l]])
    bic[l] <- -2 * loglik[l] + log(n) * Kpar
  }
  
  return(bic)
    
}