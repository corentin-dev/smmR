#' Loglikelihood
#'
#' @description Generic function computing the loglikelihood of the model `x`,
#'   with the list of sequences `sequences`.
#'
#' @param x An object for which the logikelihood can be computed.
#' @param sequences A list of vectors representing the sequences for which the 
#'   log-likelihood must be computed.
#' @param states Vector of state space (of length s).
#' @return A vector giving the value of the loglikelihood for each sequence.
#' 
#' 
#' @export
#'
loglik <- function(x, sequences, states) {
  UseMethod("loglik", x)
}