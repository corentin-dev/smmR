#' Loglikelihood
#'
#' @description Generic function computing the loglikelihood of the model `x`,
#'   with the list of sequences `seq`.
#'
#' @param x An object for which the logikelihood can be computed.
#' @param seq A list of vectors representing the sequences for which the 
#'   log-likelihood must be computed.
#' @param E Vector of state space (of length S).
#' @return A vector giving the value of the loglikelihood for each sequence.
#' 
#' 
#' @export
#'
loglik <- function(x, seq, E) {
  UseMethod("loglik", x)
}