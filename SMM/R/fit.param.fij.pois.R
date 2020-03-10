.fit.param.fij.pois <- function(res, i, j, Kmax, cens.beg) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  theta0 <- sum(0:(Kmax - 1) * res$Nijk[i, j, ]) / sum(res$Nijk[i, j, ])
  
  if (cens.beg) {# Censoring at the beginning
    
    loglik <- function(par) {
      
      mask <- res$Nijk[i, j, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      fk <- rep.int(x = 0, times = Kmax)
      fk[mask] <- dpois(x = kmask, lambda = par, log = TRUE)
      
      mask <- res$Nbijk[i, j, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = Kmax)
      Fk[mask] <- ppois(q = kmask, lambda = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(res$Nijk[i, j, ] * fk) + sum(res$Nbijk[i, j, ] * Fk)))
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = Kmax - 1)
    theta <- CO2$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(c(theta, NA))
  
}