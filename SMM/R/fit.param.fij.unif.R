.fit.param.fij.unif <- function(res, i, j, Kmax) {
  
  theta <- tail(which(res$Nijk[i, j, ] != 0), 1)
  
  return(c(theta, NA))
  
}