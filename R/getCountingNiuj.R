#' Function giving the value of the counting process Niuj used in the 
#' estimation of the kernel and the transition matrix of censored and 
#' non-parametric semi-markov chains (cf. article Exact MLE and asymptotic 
#' properties for nonparametric semi-Markov models)
#' 
#' @param sequences A list of sequences of states.
#' @return An array giving the values of the counting process \eqn{N_{iuj}}.
#' 
#' @export
#' 
getNiuj <- function(sequences) {
 
  #############################
  # Checking parameters sequences
  #############################
  
  if (!(is.list(sequences) & all(sapply(sequences, class) %in% c("character", "numeric")))) {
    stop("'sequences' should be a list of vectors")
  }
  
  states <- unique(unlist(sequences))
  processes <- processesSemiMarkov(sequences = sequences, states = states, verbose = FALSE)
  
  Ym <- lapply(processes$Ym, function(x) x - 1)
  Um <- processes$Um
  s <- processes$s
  kmax <- processes$kmax
  
  Niuj <- getCountingNiuj(Ym, Um, s, kmax)
  
  return(Niuj)
   
}
