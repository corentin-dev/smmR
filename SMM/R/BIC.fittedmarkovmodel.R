BIC.fittedmarkovmodel <- function(object, ...) {
  
  S <- length(object$estimate$E)
  
  # Kpar: number of parameters of the model
  Kpar <- (S - 1) * S ^ object$estimate$k
  
  n <- sapply(object$seq, length)
  
  vecBIC <- -2 * object$logliks + log(n) * Kpar
  
  return(vecBIC)
}