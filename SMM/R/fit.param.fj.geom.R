.fit.param.fj.geom <- function(counting, j, kmax, cens.beg) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  theta0 <- sum(counting$Njk[j, ]) / sum(1:kmax * counting$Njk[j, ])
  
  if (cens.beg) {# Censoring at the beginning
    
    loglik <- function(par) {
      
      mask <- counting$Njk[j, ] != 0
      kmask <- (0:(kmax - 1))[mask]
      fk <- rep.int(x = 0, times = kmax)
      fk[mask] <- dgeom(x = kmask, prob = par, log = TRUE)
      
      mask <- counting$Nbjk[j, ] != 0
      kmask <- (0:(kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = kmax)
      Fk[mask] <- pgeom(q = kmask, prob = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(counting$Njk[j, ] * fk) + sum(counting$Nbjk[j, ] * Fk)))
    }
    
    mle <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = 1)
    theta <- mle$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(c(theta, NA))
  
}
