.fit.param.fij.unif <- function(counting, i, j, kmax) {
  
  theta <- tail(which(counting$Nijk[i, j, ] != 0), 1)
  
  return(c(theta, NA))
  
}
