#' Loglikelihood
#'
#' @description Generic function computing the loglikelihood of the model `x`,
#'   with the list of sequences `sequences`.
#'
#' @param x An object for which the logikelihood can be computed.
#' @param sequences A list of vectors representing the sequences for which the 
#'   log-likelihood must be computed.
#' @return Value of the loglikelihood for each sequence.
#' 
#' 
#' @export
#'
loglik <- function(x, sequences) {
  UseMethod("loglik", x)
}