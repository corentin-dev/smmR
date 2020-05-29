#' Compute processes and counting processes
#'
#' @description Compute processes such as Y, T, U... and counting processes.
#' 
#' @param sequences A list of vectors representing the sequences.
#' @param states Vector of state space (of length s).
#' @return An object of S3 class processes.
#'
#' @noRd
#' 
processes <- function(sequences, states) {
  
  #############################
  # Checking parameters
  #############################
  
  if (!is.list(sequences)) {
    stop("The parameter sequences should be a list")
  }
  
  if (!all(unique(unlist(sequences)) %in% states)) {
    stop("Some states in the list of observed sequences sequences are not in the state space states")
  }
  
  s <- length(states) # State space size
  L <- length(sequences) # Number of sequences
  
  #############################
  # Get the processes
  #############################
  
  kmax <- 0 # max(max_l(n_l))
  
  Ym <- list() # Semi-Markov chain
  Jm <- list() # Successively visited states (coded with numbers)
  Tm <- list() # Successive time points when state changes
  Lm <- list() # Sojourn time
  Um <- list() # Backward recurrence time
  
  for (l in 1:L) {
    
    processes <- .getProcesses(sequences[[l]], states) # Compute the processes
    
    # Allocate the processes
    Ym[[l]] <- processes$Ym
    Jm[[l]] <- processes$Jm
    Tm[[l]] <- processes$Tm
    Lm[[l]] <- processes$Lm
    Um[[l]] <- processes$Um
    kmax <- max(kmax, max(Lm[[l]]))
  }
  
  #############################
  # Get the counting processes
  #############################
  
  counting <- .getCountingProcesses(Jm, Lm, s, kmax)
  
  
  #############################
  # Create the object sequences
  #############################
  
  ans <-
    list(
      states = states,
      s = s,
      L = L,
      Ym = Ym,
      Jm = Jm,
      Tm = Tm,
      Lm = Lm,
      Um = Um,
      kmax = kmax,
      counting = counting
    )
  
  class(ans) <- "processes"
  
  return(ans)
}

# Function to check if an object is of class processes
is.processes <- function(x) {
  inherits(x, "processes")
}