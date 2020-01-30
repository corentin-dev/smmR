fitmarkovmodel <- function(seq, E, k = 1) {
  
  #############################
  # Checking parameters seq and E
  #############################
  
  if (!(is.list(seq))) {
    stop("The parameter seq should be a list")
  }
  
  if (!all(unique(unlist(seq)) %in% E)) {
    stop("Some states in the list of observed sequences seq are not in the state space E")
  }
  
  #############################
  # Checking parameter k
  #############################
  
  if (!((k > 0) && ((k %% 1) == 0))) {
    stop("k must be a strictly positive integer")
  }
  
  
  S <- length(E) # State space size
  nbseq <- length(seq) # Number of sequences
  vect.seq <- c()
  
  # Count the number of transitions from i to j in k steps
  Nijl <- array(0, c(S ^ k, S, nbseq))
  
  for (i in 1:nbseq) {
    
    vect.seq <- seq[[i]]
    Nijl[, , i] <- matrix(count(seq = vect.seq, wordsize = k + 1, alphabet = E), byrow = TRUE, ncol = S)
    
  }
  
  Nij <- apply(Nijl, c(1, 2), sum)
  Ni <- apply(Nijl, 1, sum)
  
  # Verify the existence of all transitions
  if (length(which(Nij == 0)) > 0) {
    warning("Missing transitions")
  }
  
  # Verify the existence of all states
  if (length(which(Ni == 0)) > 0) {
    warning("Missing observed states")
  }
  
  # Compute the transition matrix
  ptrans <- Nij / tcrossprod(Ni, rep.int(1, S))
  ptrans[which(is.na(ptrans))] <- 0
  
  if (k == 1) {
    init <- .stationary.distribution(ptrans)
  } else {
    Nstart <- as.vector(count(seq = unlist(seq), wordsize = 1, alphabet = E)) ## initial law for k > 1 ???
    init <- Nstart / sum(Nstart)
  }
  
  #############################
  # Compute the log-likelihood
  #############################
  logliks <- rep.int(NA, nbseq)
  for (j in 1:nbseq) {
    s <- 0
    for (i in 1:k) {# Warning to initial law
      s <- s + log(init[which(E == seq[[j]][i])])
    }
    logliks[j] <- s + sum(as.numeric(Nijl[, , j])[which(ptrans != 0)] * log(ptrans[which(ptrans != 0)]))
  }
  
  
  estimate <- list(E = E, init = init, ptrans = ptrans, k = k)
  class(estimate) <- "markovmodel"
  
  ret <- list(estimate = estimate, seq = seq, logliks = logliks)
  class(ret) <- "fittedmarkovmodel"
  
  return(ret)
}
