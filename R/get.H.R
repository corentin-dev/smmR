# Function to compute the value of \eqn{H}
.get.H <- function(q) {
  
  k <- dim(q)[3] - 1
  
  hik <- apply(X = q, MARGIN = c(1, 3), sum)
  Hik <- t(apply(X = hik, MARGIN = 1, cumsum))
  
  H <- array(data = 0, dim = c(nrow(q), ncol(q), k + 1))
  
  for (j in 1:(k + 1)) {
    H[, , j] <- diag(Hik[, j])
  }
  
  return(H)
}


#' Function to compute the value of the sojourn time cumulative distribution \eqn{H}
#' 
#' @description Function to compute the value of \eqn{H} (See equation (3.4) p.46).
#' 
#' @param q An array giving the values of the kernel for a giving time horizon 
#'   \eqn{[0, \dots, k]} (This kernel `q` is the output of the method `getKernel` 
#'   or `get.qy`).
#' @return An array giving the value of \eqn{H(k)} at each time between 0 
#'   and `k`.
#' 
#' @export
#' 
get.H <- function(q) {
  
  #############################
  # checking parameters q
  #############################
  .is.kernel(q)
  
  H <- .get.H(q)
  
  return(H)
}
