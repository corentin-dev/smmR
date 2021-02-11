#' Log-likelihood Function
#'
#' @description Computation of the log-likelihood for a semi-Markov model
#'
#' @param x An object for which the log-likelihood can be computed.
#' @param processes An object of class `processes`.
#'
#' @noRd
#'
.loglik <- function(x, processes) {
  UseMethod(".loglik", x)
}

#' Log-likelihood Function
#'
#' @description Generic function computing the log-likelihood of the model `x`,
#'   with the list of sequences `sequences`.
#'
#' @param x An object for which the log-likelihood can be computed.
#' @param sequences A list of vectors representing the sequences for which the 
#'   log-likelihood will be computed based on `x`.
#' @return Value of the log-likelihood.
#' 
#' @export
#'
loglik <- function(x, sequences) {
  UseMethod("loglik", x)
}
