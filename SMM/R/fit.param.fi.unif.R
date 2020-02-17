.fit.param.fi.unif <- function(res, i, Kmax, cens.beg, cens.end) {
  
  theta <- tail(which(res$Nik[i, ] != 0), 1)
  
  return(c(theta, NA))
  
}