.fit.param.fj.pois <- function(res, j, Kmax, cens.beg) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  theta0 <- sum(0:(Kmax - 1) * res$Njk[j, ]) / sum(res$Njk[j, ])
  
  if (cens.beg) {# Censoring at the beginning
    
    loglik <- function(par) {
      
      mask <- res$Njk[j, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      fk <- rep.int(x = 0, times = Kmax)
      fk[mask] <- dpois(x = kmask, lambda = par, log = TRUE)
      
      mask <- res$Nbjk[j, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = Kmax)
      Fk[mask] <- ppois(q = kmask, lambda = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(res$Njk[j, ] * fk) + sum(res$Nbjk[j, ] * Fk)))
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = Kmax - 1)
    theta <- CO2$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(c(theta, NA))
  
}