.fit.param.fj.pois <- function(counting, j, kmax, cens.beg) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  theta0 <- sum(0:(kmax - 1) * counting$Njk[j, ]) / sum(counting$Njk[j, ])
  
  if (cens.beg) {# Censoring at the beginning
    
    loglik <- function(par) {
      
      mask <- counting$Njk[j, ] != 0
      kmask <- (0:(kmax - 1))[mask]
      fk <- rep.int(x = 0, times = kmax)
      fk[mask] <- dpois(x = kmask, lambda = par, log = TRUE)
      
      mask <- counting$Nbjk[j, ] != 0
      kmask <- (0:(kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = kmax)
      Fk[mask] <- ppois(q = kmask, lambda = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(counting$Njk[j, ] * fk) + sum(counting$Nbjk[j, ] * Fk)))
    }
    
    mle <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = kmax - 1)
    theta <- mle$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(c(theta, NA))
  
}
