.fit.param.f.geom <- function(counting, kmax, cens.beg, cens.end) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  theta0 <- sum(counting$Nk) / sum(1:kmax * counting$Nk)
  
  if (!cens.beg & cens.end) {# Censoring at the end
    
    loglik <- function(par) {
      
      mask <- counting$Nk != 0
      kmask <- (0:(kmax - 1))[mask]
      fk <- rep.int(x = 0, times = kmax)
      fk[mask] <- dgeom(x = kmask, prob = par, log = TRUE)
      
      mask <- counting$Nek != 0
      kmask <- (0:(kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = kmax)
      Fk[mask] <- pgeom(q = kmask, prob = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(counting$Nk * fk) + sum(counting$Nek * Fk)))
    }
    
    mle <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = 1)
    theta <- mle$par
    
  } else if (cens.beg & !cens.end) {# censoring at the beginning
    
    loglik <- function(par) {
      
      mask <- counting$Nk != 0
      kmask <- (0:(kmax - 1))[mask]
      fk <- rep.int(x = 0, times = kmax)
      fk[mask] <- dgeom(x = kmask, prob = par, log = TRUE)
      
      mask <- counting$Nbk != 0
      kmask <- (0:(kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = kmax)
      Fk[mask] <- pgeom(q = kmask, prob = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(counting$Nk * fk) + sum(counting$Nbk * Fk)))
    }
    
    mle <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = 1)
    theta <- mle$par
    
  } else if (cens.beg & cens.end) {# censoring at the beginningand at the end
    
    loglik <- function(par) {
      
      mask <- counting$Nk != 0
      kmask <- (0:(kmax - 1))[mask]
      fk <- rep.int(x = 0, times = kmax)
      fk[mask] <- dgeom(x = kmask, prob = par, log = TRUE)
      
      mask <- counting$Nebk != 0
      kmask <- (0:(kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = kmax)
      Fk[mask] <- pgeom(q = kmask, prob = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(counting$Nk * fk) + sum(counting$Nebk * Fk)))
    }
    
    mle <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = 1)
    theta <- mle$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(c(theta, NA))
  
}
