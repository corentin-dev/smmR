# Function to compute the value of \eqn{\psi}
.get.psi <- function(q) {
  
  k <- dim(q)[3] - 1
  
  psi <- array(data = 0, dim = c(nrow(q), ncol(q), k + 1)) # (S, S, k + 1)
  psi[, , 1] <- diag(x = 1, nrow = nrow(q), ncol = ncol(q)) # k = 0
  
  for (j in 1:k) {
    
    psi[, , j + 1] <-
      -Reduce('+', lapply(
        X = 0:(j - 1),
        FUN = function(l)
          psi[, , l + 1] %*% (-q[, , j - l + 1])
      ))
  }
  
  return(psi)
  
}

#' Function to compute the value of the matrix-valued function \eqn{\psi}
#' 
#' @description Function to compute the value of \eqn{\psi}, the matrix-valued 
#'   function (See equation (3.16) p.53).
#' 
#' @param q An array giving the values of the kernel for a giving time horizon 
#'   \eqn{[0, \dots, k]} (This kernel `q` is the output of the method `getKernel` 
#'   or `.get.qy`).
#' @return An array giving the value of \eqn{\psi(k)} at each time between 0 
#'   and `k`.
#' 
#' @export
#' 
get.psi <- function(q) {
  
  .is.kernel(q)
  
  psi <- .get.psi(q)
  
  return(psi)
  
}
