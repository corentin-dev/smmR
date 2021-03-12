#' Compute processes and counting processes for semi-Markov models
#' 
#' @description Compute processes such as Y, T, U... and counting processes.
#' 
#' @param sequences A list of vectors representing the sequences.
#' @param states Vector of state space (of length s).
#' @param verbose Boolean. If TRUE, messages are generated.
#' @return An object of class S3 `processesSemiMarkov`.
#' 
#' @noRd
#' 
processesSemiMarkov <- function(sequences, states, verbose = TRUE) {
  
  s <- length(states) # State space size
  L <- length(sequences) # Number of sequences
  
  #############################
  # Get the processes
  #############################
  
  numericSequences <- lapply(sequences, function(x) sapply(x, function(y) which(states == y) - 1, USE.NAMES = FALSE))
  processes <- getProcesses(numericSequences) # Compute the processes
  
  Ym <- lapply(processes$Ym, function(x) x + 1) # Semi-Markov chain
  Jm <- lapply(processes$Jm, function(x) x + 1) # Successively visited states (coded with numbers)
  Tm <- processes$Tm # Successive time points when state changes
  Lm <- processes$Lm # Sojourn time
  Um <- processes$Um # Backward recurrence time
  M <- sum(sapply(processes$Ym, length))
  kmax <- max(unlist(Lm))
  
  #############################
  # Get the counting processes
  #############################
  
  counting <- getCountingProcesses(lapply(Jm, function(x) x - 1), Lm, s, kmax)
  
  if (verbose) {
    indexdiag <- seq(1, s * s, by = s + 1)
    Nij <- as.vector(counting$Nij)[-indexdiag]
    statesi <- .row(c(s, s))[-indexdiag][which(Nij == 0)]
    statesj <- .col(c(s, s))[-indexdiag][which(Nij == 0)]
    
    if ((length(statesi) != 0) | (length(statesj) != 0)) {
      message(
        "Some transitions from state i to state j are not observed: ",
        paste0(sapply(1:length(statesi), function(x)
          paste0("(i = \"", states[statesi[x]], "\" to j = \"", states[statesj[x]], "\")")), collapse = ", ")
      )
    }
  }
  
  #############################
  # Create the object sequences
  #############################
  
  ans <-
    list(
      states = states,
      s = s,
      M = M,
      L = L,
      kmax = kmax,
      Ym = Ym,
      Jm = Jm,
      Tm = Tm,
      Lm = Lm,
      Um = Um,
      counting = counting
    )
  
  class(ans) <- "processesSemiMarkov"
  
  return(ans)
}


#' Compute processes and counting processes for Markov models
#' 
#' @description Compute counting processes.
#' 
#' @param sequences A list of vectors representing the sequences.
#' @param states Vector of state space (of length s).
#' @param verbose Boolean. If TRUE, messages are generated.
#' @return An object of class S3 `processesMarkov`.
#' 
#' @noRd
#' 
processesMarkov <- function(sequences, states, k, verbose = TRUE) {
  
  s <- length(states) # State space size
  L <- length(sequences) # Number of sequences
  M <- sum(sapply(sequences, length))
  
  # Count the number of transitions from i to j in k steps
  Nijl <- array(0, c(s ^ k, s, L))
  Nstarti <- rep.int(x = 0, times = s ^ k)
  
  for (i in 1:L) {
    
    Nijl[, , i] <- matrix(count(seq = sequences[[i]], wordsize = k + 1, alphabet = states), byrow = TRUE, ncol = s)
    
    Nstarti[which(words(length = k, alphabet = states) == c2s(sequences[[i]][1:k]))] <- Nstarti[which(words(length = k, alphabet = states) == c2s(sequences[[i]][1:k]))] + 1
    
  }
  
  Nij <- apply(Nijl, c(1, 2), sum)
  Ni <- apply(Nijl, 1, sum)
  
  if (verbose) {
    
    # Verify the existence of all states
    if (length(which(Ni == 0)) > 0) {
      warning("Missing observed states: \"", paste0(states[which(Ni == 0)]), "\"")
    }
    
  }
  
  
  #############################
  # Create the object sequences
  #############################
  
  ans <-
    list(
      states = states,
      s = s,
      L = L,
      M = M, 
      Nij = Nij, 
      Ni = Ni, 
      Nstarti = Nstarti
    )
  
  class(ans) <- "processesMarkov"
  
  return(ans)
}
