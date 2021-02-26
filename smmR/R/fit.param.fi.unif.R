.fit.param.fi.unif <- function(counting, i, kmax) {
  
  theta <- tail(which(counting$Nik[i, ] != 0), 1)
  
  return(c(theta, NA))
  
}
