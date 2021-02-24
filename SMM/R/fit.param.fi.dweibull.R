.fit.param.fi.dweibull <- function(counting, i, kmax, cens.beg, cens.end) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  theta0 <- suppressWarnings(estdweibull(x = unlist(sapply(1:kmax, function(x) rep(x, counting$Nik[i, x]))), method = "ML", zero = FALSE))
  
  loglik <- function(par) {
    
    mask <- counting$Nik[i, ] != 0
    kmask <- (1:kmax)[mask]
    fk <- rep.int(x = 0, times = kmax)
    fk[mask] <- log(ddweibull(x = kmask, q = par[1], beta = par[2], zero = FALSE))
    
    return(-(sum(counting$Nik[i, ] * fk)))
  }
  
  # Constraints about the values of the parameters:
  
  # q, beta > 0
  u0 <- diag(x = 1, nrow = 2)
  c0 <- c(0, 0)
  
  # q < 1
  u1 <- matrix(data = c(-1, 0), nrow = 1, ncol = 2)
  c1 <- c(-1)
  
  mle <- constrOptim(
    theta = theta0,
    f = loglik,
    ui = rbind(u0, u1),
    ci = c(c0, c1),
    method = "Nelder-Mead"
  )
  theta0 <- mle$par
  
  if (!cens.beg & cens.end) {# Censoring at the end
    
    loglik <- function(par) {
      
      mask <- counting$Nik[i, ] != 0
      kmask <- (1:kmax)[mask]
      fk <- rep.int(x = 0, times = kmax)
      fk[mask] <- log(ddweibull(x = kmask, q = par[1], beta = par[2], zero = FALSE))
      
      mask <- counting$Neik[i, ] != 0
      kmask <- (1:kmax)[mask]
      Fk <- rep.int(x = 0, times = kmax)
      Fk[mask] <- log(1 - pdweibull(x = kmask, q = par[1], beta = par[2], zero = FALSE))
      
      return(-(sum(counting$Nik[i, ] * fk) + sum(counting$Neik[i, ] * Fk)))
    }
    
    mle <- constrOptim(
      theta = theta0,
      f = loglik,
      ui = rbind(u0, u1),
      ci = c(c0, c1),
      method = "Nelder-Mead"
    )
    theta <- mle$par
    
  } else if (cens.beg & !cens.end) {# Censoring at the beginning
    
    loglik <- function(par) {
      
      mask <- counting$Nik[i, ] != 0
      kmask <- (1:kmax)[mask]
      fk <- rep.int(x = 0, times = kmax)
      fk[mask] <- log(ddweibull(x = kmask, q = par[1], beta = par[2], zero = FALSE))
      
      mask <- counting$Nbik[i, ] != 0
      kmask <- (1:kmax)[mask]
      Fk <- rep.int(x = 0, times = kmax)
      Fk[mask] <- log(1 - pdweibull(x = kmask, q = par[1], beta = par[2], zero = FALSE))
      
      return(-(sum(counting$Nik[i, ] * fk) + sum(counting$Nbik[i, ] * Fk)))
    }
    
    mle <- constrOptim(
      theta = theta0,
      f = loglik,
      ui = rbind(u0, u1),
      ci = c(c0, c1),
      method = "Nelder-Mead"
    )
    theta <- mle$par
    
  } else if (cens.beg & cens.end) {# Censoring at the beginning and at the end
    
    loglik <- function(par) {
      
      mask <- counting$Nik[i, ] != 0
      kmask <- (1:kmax)[mask]
      fk <- rep.int(x = 0, times = kmax)
      fk[mask] <- log(ddweibull(x = kmask, q = par[1], beta = par[2], zero = FALSE))
      
      mask <- counting$Nebik[i, ] != 0
      kmask <- (1:kmax)[mask]
      Fk <- rep.int(x = 0, times = kmax)
      Fk[mask] <- log(1 - pdweibull(x = kmask, q = par[1], beta = par[2], zero = FALSE))
      
      return(-(sum(counting$Nik[i, ] * fk) + sum(counting$Nebik[i, ] * Fk)))
    }
    
    mle <- constrOptim(
      theta = theta0,
      f = loglik,
      ui = rbind(u0, u1),
      ci = c(c0, c1),
      method = "Nelder-Mead"
    )
    theta <- mle$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(theta)
  
}
