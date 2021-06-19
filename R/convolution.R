#' Discrete-time convolution product of \eqn{f} and \eqn{g}
#'   (See definition 2.2 p. 20)
#'
#' @param f A vector giving the values of the function \eqn{f} for each
#'   \eqn{k \in \mathbb{N}}.
#' @param g A vector giving the values of the function \eqn{g} for each
#'   \eqn{k \in \mathbb{N}}.
#'
#' @return A vector giving the values of the discrete-time convolution of
#'   \eqn{f} and \eqn{g} for each \eqn{k \in \mathbb{N}}.
#'
#' @export
#'
convolution <- function(f, g) {
  #############################
  # Checking parameter f
  #############################
  
  if (!(is.numeric(f) & !anyNA(f) & is.vector(f))) {
    stop("'f' is not a numeric vector")
  }
  
  #############################
  # Checking parameter g
  #############################
  
  if (!(is.numeric(g) & !anyNA(g) & is.vector(g))) {
    stop("'g' is not a numeric vector")
  }
  
  if (!(length(f) == length(g))) {
    stop("'f' and 'g' must have the same length")
  }
  
  conv <- C_convolution(f, g)
  
  return(conv)
}
