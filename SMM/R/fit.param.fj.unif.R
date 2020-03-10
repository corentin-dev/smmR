.fit.param.fj.unif <- function(res, j, Kmax) {
  
  theta <- tail(which(res$Njk[j, ] != 0), 1)
  
  return(c(theta, NA))
  
}