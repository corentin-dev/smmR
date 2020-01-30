AIC.fittedmarkovmodel <- function(object, ..., k = 2) {
  
  S <- length(object$estimate$E)
  
  # Kpar: number of parameters of the model
  Kpar <- (S - 1) * S ^ object$estimate$k
  
  vecAIC <- -2 * object$logliks + k * Kpar
  
  return(vecAIC)
  
}