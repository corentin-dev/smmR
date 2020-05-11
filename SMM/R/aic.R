#' Akaike Information Criterion (AIC)
#'
#' @description Generic function computing the Akaike Information Criterion of 
#'   the model `x`, with the list of sequences `seq`.
#'
#' @param x An object for which the logikelihood can be computed.
#' @param seq A list of vectors representing the sequences for which the 
#'   AIC criterion must be computed.
#' @param E Vector of state space (of length S).
#' @return A numeric value giving the value of the AIC.
#' 
#' 
#' @export
#'
aic <- function(x, seq, E) {
  UseMethod("aic", x)
}