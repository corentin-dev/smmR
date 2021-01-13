#' Compute processes and counting processes
#'
#' @description Compute processes such as Y, T, U... and counting processes.
#' 
#' @param sequences A list of vectors representing the sequences.
#' @param states Vector of state space (of length s).
#' @param verbose Boolean. If TRUE, messages are generated.
#' @return An object of class S3 `processes`.
#'
#' @noRd
#' 
processes <- function(sequences, states, verbose = TRUE) {
  
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
  
  numericSequences <- lapply(sequences, function(x) sapply(x, function(y) which(states == y) - 1, USE.NAMES = FALSE))
  processes <- getProcesses(numericSequences) # Compute the processes
  
  Ym <- lapply(processes$Ym, function(x) x + 1) # Semi-Markov chain
  Jm <- lapply(processes$Jm, function(x) x + 1) # Successively visited states (coded with numbers)
  Tm <- processes$Tm # Successive time points when state changes
  Lm <- processes$Lm # Sojourn time
  Um <- processes$Um # Backward recurrence time
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
    
    if ((length(statesi) != 0) |
        (length(statesj) != 0)) {
      message(
        "Some transitions from state i to state j are not observed: ",
        paste0(sapply(1:length(statesi), function(x)
          paste0("(i=", statesi[x], " to j=", statesj[x], ")")), collapse = ", ")
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
