#' Discrete-time matrix convolution product
#'   (See definition 3.5 p. 48)
#'
#' @param a An array of dimension \eqn{(S, S, k + 1)}.
#' @param b An array of dimension \eqn{(S, S, k + 1)}.
#'
#' @return An array of dimension \eqn{(S, S, k + 1)} giving the discrete-time
#'   matrix convolution product for each \eqn{k \in \mathbb{N}}.
#'
#' @export
#'
matrixConvolution <- function(a, b) {
  #############################
  # Checking parameter a
  #############################
  
  if (!(is.numeric(a) & !anyNA(a) & is.array(a) & !is.matrix(a) & (dim(a)[1] == dim(a)[2]))) {
    stop("'a' must be a numeric array of dimension (s, s, k+1)")
  }
  
  #############################
  # Checking parameter b
  #############################
  
  if (!(is.numeric(b) & !anyNA(b) & is.array(b) & !is.matrix(b) & (dim(b)[1] == dim(b)[2]))) {
    stop("'b' must be a numeric array of dimension (s, s, k+1)")
  }
  
  if (!((dim(a)[1] == dim(b)[1]) & (dim(a)[3] == dim(b)[3]))) {
    stop("'a' and 'b' must have the same dimensions")
  }
  
  c <- C_matrixConvolution(a, b)
  
  return(c)
}
