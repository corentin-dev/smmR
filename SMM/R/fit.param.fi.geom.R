.fit.param.fi.geom <- function(res, i, Kmax, cens.beg, cens.end) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  theta0 <- sum(res$Nik[i, ]) / sum(1:Kmax * res$Nik[i, ])

  if (!cens.beg && cens.end) {# Censoring at the end
    
    loglik <- function(par) {
      
      mask <- res$Nik[i, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      fk <- rep.int(x = 0, times = Kmax)
      fk[mask] <- dgeom(x = kmask, prob = par, log = TRUE)
      
      mask <- res$Neik[i, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = Kmax)
      Fk[mask] <- pgeom(q = kmask, prob = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(res$Nik[i, ] * fk) + sum(res$Neik[i, ] * Fk)))
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = 1)
    theta <- CO2$par

  } else if (cens.beg && !cens.end) {# censoring at the beginning
    
    loglik <- function(par) {
      
      mask <- res$Nik[i, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      fk <- rep.int(x = 0, times = Kmax)
      fk[mask] <- dgeom(x = kmask, prob = par, log = TRUE)
      
      mask <- res$Nbik[i, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = Kmax)
      Fk[mask] <- pgeom(q = kmask, prob = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(res$Nik[i, ] * fk) + sum(res$Nbik[i, ] * Fk)))
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = 1)
    theta <- CO2$par
    
  } else if (cens.beg && cens.end) {# censoring at the beginningand at the end
    
    loglik <- function(par) {
      
      mask <- res$Nik[i, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      fk <- rep.int(x = 0, times = Kmax)
      fk[mask] <- dgeom(x = kmask, prob = par, log = TRUE)
      
      mask <- res$Nebik[i, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = Kmax)
      Fk[mask] <- pgeom(q = kmask, prob = par, lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(res$Nik[i, ] * fk) + sum(res$Nebik[i, ] * Fk)))
    }

    CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = 1)
    theta <- CO2$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(c(theta, NA))
    
}