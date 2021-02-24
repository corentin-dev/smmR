.fit.param.fi.pois <- function(counting, i, kmax, cens.beg, cens.end) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  theta0 <- sum(0:(kmax - 1) * counting$Nik[i, ]) / sum(counting$Nik[i, ])
  
  if (!cens.beg & cens.end) {# Censoring at the end
    
    loglik <- function(par) {
      
      mask <- counting$Nik[i, ] != 0
      kmask <- (0:(kmax - 1))[mask]
      fk <- rep.int(x = 0, times = kmax)
      fk[mask] <- dpois(x = kmask, lambda = par, log = TRUE)
      
      mask <- counting$Neik[i, ] != 0
      kmask <- (0:(kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = kmax)
      Fk[mask] <- ppois(q = kmask, lambda = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(counting$Nik[i, ] * fk) + sum(counting$Neik[i, ] * Fk)))
    }
    
    mle <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = kmax - 1)
    theta <- mle$par
    
  } else if (cens.beg & !cens.end) {# Censoring at the beginning
    
    loglik <- function(par) {
      
      mask <- counting$Nik[i, ] != 0
      kmask <- (0:(kmax - 1))[mask]
      fk <- rep.int(x = 0, times = kmax)
      fk[mask] <- dpois(x = kmask, lambda = par, log = TRUE)
      
      mask <- counting$Nbik[i, ] != 0
      kmask <- (0:(kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = kmax)
      Fk[mask] <- ppois(q = kmask, lambda = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(counting$Nik[i, ] * fk) + sum(counting$Nbik[i, ] * Fk)))
    }
    
    mle <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = kmax - 1)
    theta <- mle$par
    
  } else if (cens.beg & cens.end) {# Censoring at the beginning and at the end
    
    loglik <- function(par) {
      
      mask <- counting$Nik[i, ] != 0
      kmask <- (0:(kmax - 1))[mask]
      fk <- rep.int(x = 0, times = kmax)
      fk[mask] <- dpois(x = kmask, lambda = par, log = TRUE)
      
      mask <- counting$Nebik[i, ] != 0
      kmask <- (0:(kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = kmax)
      Fk[mask] <- ppois(q = kmask, lambda = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(counting$Nik[i, ] * fk) + sum(counting$Nebik[i, ] * Fk)))
    }
    
    mle <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = kmax - 1)
    theta <- mle$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(c(theta, NA))
  
}
