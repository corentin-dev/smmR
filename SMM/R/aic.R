#' Akaike Information Criterion (AIC)
#'
#' @description Generic function computing the Akaike Information Criterion of 
#'   the model `x`, with the list of sequences `sequences`.
#'
#' @param x An object for which the log-likelihood can be computed.
#' @param sequences A list of vectors representing the sequences for which the 
#'   AIC will be computed based on `x`.
#' @return Value of the AIC.
#' 
#' @export
#'
aic <- function(x, sequences) {
  UseMethod("aic", x)
}
