#' Compute processes and counting processes
#'
#' @description Compute processes such as Y, T, U... and counting processes.
#' 
#' @param seq A list of vectors representing the sequences.
#' @param E Vector of state space (of length S).
#' @return An object of S3 class sequences.
#'
#' @noRd
#' 
sequences <- function(seq, E) {
  
  #############################
  # Checking parameters
  #############################
  
  if (!is.list(seq)) {
    stop("The parameter seq should be a list")
  }
  
  if (!all(unique(unlist(seq)) %in% E)) {
    stop("Some states in the list of observed sequences seq are not in the state space E")
  }
  
  S <- length(E) # State space size
  L <- length(seq) # Number of sequences
  
  #############################
  # Get the processes
  #############################
  
  Kmax <- 0 # max(max_l(n_l))
  
  Ym <- list() # Semi-Markov chain
  Jm <- list() # Successively visited states (coded with numbers)
  Tm <- list() # Successive time points when state changes
  Lm <- list() # Sojourn time
  Um <- list() # Backward recurrence time
  
  for (l in 1:L) {
    
    processes <- .getProcesses(seq[[l]], E) # Compute the processes
    
    # Allocate the processes
    Ym[[l]] <- processes$Ym
    Jm[[l]] <- processes$Jm
    Tm[[l]] <- processes$Tm
    Lm[[l]] <- processes$Lm
    Um[[l]] <- processes$Um
    Kmax <- max(Kmax, max(Lm[[l]]))
  }
  
  #############################
  # Get the counting processes
  #############################
  
  counting <- .getCountingProcesses(Jm, Lm, S, Kmax)
  
  
  #############################
  # Create the object sequences
  #############################
  
  ans <-
    list(
      E = E,
      S = S,
      L = L,
      Ym = Ym,
      Jm = Jm,
      Tm = Tm,
      Lm = Lm,
      Um = Um,
      Kmax = Kmax,
      counting = counting
    )
  
  class(ans) <- "sequences"
  
  return(ans)
}

# Function to check if an object is of class sequences
is.sequences <- function(x) {
  inherits(x, "sequences")
}