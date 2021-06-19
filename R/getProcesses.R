#' Function to compute processes based on a list of sequences
#' 
#' @param sequences A list of sequences of states.
#' @return A list giving:
#'  \itemize{
#'    \item{s : }{Cardinal of the states space;}
#'    \item{states : }{Vector of state space of length \eqn{s};}
#'    \item{M : }{The total length (sum of the length of each sequence);}
#'    \item{L : }{The number of sequences;}
#'    \item{kmax : }{Maximum length of the sojourn times;}
#'    \item{Ym: }{List of semi-Markov chains;}
#'    \item{Jm: }{List of successively visited states for each sequence;}
#'    \item{Tm: }{List of successive time points when state changes for each sequence;}
#'    \item{Lm: }{List of sojourn time for each sequence;}
#'    \item{Um: }{List of backward recurrence time for each sequence;}
#'    \item{counting : }{List giving the counting processes such as 
#'      \eqn{N_{ij}(k)}, \eqn{N_{i}(k)}, \eqn{N_{j}(k)}...}
#'  }
#' 
#' @export
#' 
getProcesses <- function(sequences) {
  
  #############################
  # Checking parameters sequences
  #############################
  
  if (!(is.list(sequences) & all(sapply(sequences, class) %in% c("character", "numeric")))) {
    stop("'sequences' should be a list of vectors")
  }
  
  
  states <- unique(unlist(sequences))
  processes <- processesSemiMarkov(sequences = sequences, states = states, verbose = FALSE)
  
  processes$Ym <- sequences
  processes$Jm <- lapply(processes$Jm, function(x) states[x])
  class(processes) <- NULL
  
  return(processes)
  
}
