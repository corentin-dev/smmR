sequences <- function(listSeq, E) {
  
  #############################
  # Checking parameters
  #############################
  
  if (!is.list(listSeq)) {
    stop("The parameter listSeq should be a list")
  }
  
  if (!all(unique(unlist(listSeq)) %in% E)) {
    stop("Some states in the list of observed sequences listSeq are not in the state space E")
  }
  
  S <- length(E)
  nbSeq <- length(listSeq)
  
  #############################
  # Get the processes
  #############################
  
  Kmax <- 0
  
  Y <- list() # Semi-Markov chain
  Jchar <- listSeq # Successively visited states
  J <- list() # Successively visited states (coded with numbers)
  T <- list() # Successive time points when state changes
  L <- list() # Sojourn time
  U <- list() # Backward recurrence time
  
  for (m in 1:nbSeq) {
    
    processes <- .getProcesses(listSeq[[m]], E)
    
    Y[[m]] <- processes$Y
    J[[m]] <- processes$J
    T[[m]] <- processes$T
    L[[m]] <- processes$L
    U[[m]] <- processes$U
    Kmax <- max(Kmax, max(L[[m]]))
  }
  
  #############################
  # Get the counting processes
  #############################
  
  counting <- .getCountingProcesses(J, L, S, Kmax)
  
  
  ans <-
    list(
      E = E,
      S = S,
      nbSeq = nbSeq,
      Y = Y,
      J = J,
      T = T,
      L = L,
      U = U,
      Kmax = Kmax,
      counting = counting
    )
  
  class(ans) <- "sequences"
  
  return(ans)
}

is.sequences <- function(x) {
  inherits(x, "sequences")
}