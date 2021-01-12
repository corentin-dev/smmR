#' Akaike Information Criterion (AIC)
#'
#' @description Generic function computing the Akaike Information Criterion of 
#'   the model `x`, with the list of sequences `sequences`.
#'
#' @param x An object for which the logikelihood can be computed.
#' @param sequences A list of vectors representing the sequences for which the 
#'   AIC criterion must be computed.
#' @return Value of the AIC.
#' 
#' 
#' @export
#'
aic <- function(x, sequences) {
  UseMethod("aic", x)
}