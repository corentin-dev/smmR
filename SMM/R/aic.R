#' Akaike Information Criterion (AIC)
#'
#' @description Computation of the Akaike Information Criterion.
#'
#' @param x An object of class [markovmodel] or inheriting from the class smm.
#'   ([smmnonparametric][smmnonparametric] or [smmparametric][smmparametric]).
#' @param seq A list of vectors representing the sequences for which the 
#'   AIC criterion must be computed.
#' @param E Vector of state space (of length S).
#' @return A vector giving the value of the AIC for each sequence.
#' 
#' 
#' @export
#'
aic <- function(x, seq, E) {
  UseMethod("aic", x)
}