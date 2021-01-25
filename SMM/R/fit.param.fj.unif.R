.fit.param.fj.unif <- function(counting, j, kmax) {
  
  theta <- tail(which(counting$Njk[j, ] != 0), 1)
  
  return(c(theta, NA))
  
}
