#' Bayesian Information Criterion (BIC)
#'
#' @description Computation of the Bayesian Information Criterion.
#'
#' @param x An object of class [markovmodel] or inheriting from the class smm.
#'   ([smmnonparametric][smmnonparametric] or [smmparametric][smmparametric]).
#' @param seq A list of vectors representing the sequences for which the 
#'   BIC criterion must be computed.
#' @param E Vector of state space (of length S).
#' @return A vector giving the value of the BIC for each sequence.
#' 
#' 
#' @export
#'
bic <- function(x, seq, E) {
  UseMethod("bic", x)
}