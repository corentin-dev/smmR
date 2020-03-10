.fit.param.fi.unif <- function(res, i, Kmax) {
  
  theta <- tail(which(res$Nik[i, ] != 0), 1)
  
  return(c(theta, NA))
  
}