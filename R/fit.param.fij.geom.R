.fit.param.fij.geom <- function(counting, i, j, kmax, cens.beg) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  theta0 <- sum(counting$Nijk[i, j, ]) / sum(1:kmax * counting$Nijk[i, j, ])
  
  if (cens.beg) {# Censoring at the beginning
    
    logLik <- function(par) {
      
      mask <- counting$Nijk[i, j, ] != 0
      kmask <- (0:(kmax - 1))[mask]
      fk <- rep.int(x = 0, times = kmax)
      fk[mask] <- dgeom(x = kmask, prob = par, log = TRUE)
      
      mask <- counting$Nbijk[i, j, ] != 0
      kmask <- (0:(kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = kmax)
      Fk[mask] <- pgeom(q = kmask, prob = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(counting$Nijk[i, j, ] * fk) + sum(counting$Nbijk[i, j, ] * Fk)))
    }
    
    mle <- optim(par = theta0, logLik, method = "Brent", lower = 0, upper = 1)
    theta <- mle$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(c(theta, NA))
  
}
