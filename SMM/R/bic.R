#' Bayesian Information Criterion (BIC)
#'
#' @description Generic function computing the Bayesian Information Criterion 
#'   of the model `x`, with the list of sequences `sequences`.
#'
#' @param x An object for which the logikelihood can be computed.
#' @param sequences A list of vectors representing the sequences for which the 
#'   BIC criterion must be computed.
#' @param states Vector of state space E (of length s).
#' @return A numeric value giving the value of the BIC.
#' 
#' 
#' @export
#'
bic <- function(x, sequences, states) {
  UseMethod("bic", x)
}