loglik <- function(x, listSeq, E) {
  
  if (!(is.smm(x))) {
    stop("x is not an object of class smm")
  }
  
  out <- .loglik(x, listSeq, E)
  
  return(out$loglik)
  
}