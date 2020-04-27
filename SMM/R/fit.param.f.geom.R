.fit.param.f.geom <- function(res, Kmax, cens.beg, cens.end) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  theta0 <- sum(res$Nk) / sum(1:Kmax * res$Nk)

  if (!cens.beg && cens.end) {# Censoring at the end
    
    loglik <- function(par) {
      
      mask <- res$Nk != 0
      kmask <- (0:(Kmax - 1))[mask]
      fk <- rep.int(x = 0, times = Kmax)
      fk[mask] <- dgeom(x = kmask, prob = par, log = TRUE)
      
      mask <- res$Nek != 0
      kmask <- (0:(Kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = Kmax)
      Fk[mask] <- pgeom(q = kmask, prob = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(res$Nk * fk) + sum(res$Nek * Fk)))
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = 1)
    theta <- CO2$par

  } else if (cens.beg && !cens.end) {# censoring at the beginning
    
    loglik <- function(par) {
      
      mask <- res$Nk != 0
      kmask <- (0:(Kmax - 1))[mask]
      fk <- rep.int(x = 0, times = Kmax)
      fk[mask] <- dgeom(x = kmask, prob = par, log = TRUE)
      
      mask <- res$Nbk != 0
      kmask <- (0:(Kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = Kmax)
      Fk[mask] <- pgeom(q = kmask, prob = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(res$Nk * fk) + sum(res$Nbk * Fk)))
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = 1)
    theta <- CO2$par
    
  } else if (cens.beg && cens.end) {# censoring at the beginningand at the end
    
    loglik <- function(par) {
      
      mask <- res$Nk != 0
      kmask <- (0:(Kmax - 1))[mask]
      fk <- rep.int(x = 0, times = Kmax)
      fk[mask] <- dgeom(x = kmask, prob = par, log = TRUE)
      
      mask <- res$Nebk != 0
      kmask <- (0:(Kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = Kmax)
      Fk[mask] <- pgeom(q = kmask, prob = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(res$Nk * fk) + sum(res$Nebk * Fk)))
    }

    CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = 1)
    theta <- CO2$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(c(theta, NA))
    
}