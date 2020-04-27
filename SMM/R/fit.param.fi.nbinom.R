.fit.param.fi.nbinom <- function(res, i, Kmax, cens.beg, cens.end) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  xbar <- sum(0:(Kmax - 1) * res$Nik[i, ]) / sum(res$Nik[i, ])
  s2 <- (1 / (sum(res$Nik[i, ]) - 1)) * sum(res$Nik[i, ] * (0:(Kmax - 1) - xbar) ^ 2)
  
  if (xbar >= s2) {
    stop("The negative binomial is not the appropriate distribution for modeling the data")
  }
  
  alphahat <- xbar ^ 2 / (s2 - xbar)
  phat <- xbar / s2
  
  theta0 <- c(alphahat, phat)
  
  loglik <- function(par) {
    
    mask <- res$Nik[i, ] != 0
    kmask <- (0:(Kmax - 1))[mask]
    fk <- rep.int(x = 0, times = Kmax)
    fk[mask] <- dnbinom(x = kmask, size = par[1], prob = par[2], log = TRUE)
    
    return(-(sum(res$Nik[i, ] * fk)))
  }
  
  # Constraints about the values of the parameters:
  
  # alpha, p > 0
  u0 <- diag(x = 1, nrow = 2)
  c0 <- c(0, 0)
  
  # p < 1
  u1 <- matrix(data = c(0, -1), nrow = 1, ncol = 2)
  c1 <- c(-1)
  
  CO2 <- constrOptim(
    theta = theta0,
    f = loglik,
    ui = rbind(u0, u1),
    ci = c(c0, c1),
    method = "Nelder-Mead"
  )
  theta0 <- CO2$par
  
  if (!cens.beg && cens.end) {# Censoring at the end
    
    loglik <- function(par) {
      
      mask <- res$Nik[i, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      fk <- rep.int(x = 0, times = Kmax)
      fk[mask] <- dnbinom(x = kmask, size = par[1], prob = par[2], log = TRUE)
      
      mask <- res$Neik[i, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = Kmax)
      Fk[mask] <- pnbinom(q = kmask, size = par[1], prob = par[2], lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(res$Nik[i, ] * fk) + sum(res$Neik[i, ] * Fk)))
    }
    
    CO2 <- constrOptim(
      theta = theta0,
      f = loglik,
      ui = rbind(u0, u1),
      ci = c(c0, c1),
      method = "Nelder-Mead"
    )
    theta <- CO2$par
    
  } else if (cens.beg && !cens.end) {# Censoring at the beginning
    
    loglik <- function(par) {
      
      mask <- res$Nik[i, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      fk <- rep.int(x = 0, times = Kmax)
      fk[mask] <- dnbinom(x = kmask, size = par[1], prob = par[2], log = TRUE)
      
      mask <- res$Nbik[i, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = Kmax)
      Fk[mask] <- pnbinom(q = kmask, size = par[1], prob = par[2], lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(res$Nik[i, ] * fk) + sum(res$Nbik[i, ] * Fk)))
    }
    
    CO2 <- constrOptim(
      theta = theta0,
      f = loglik,
      ui = rbind(u0, u1),
      ci = c(c0, c1),
      method = "Nelder-Mead"
    )
    theta <- CO2$par
    
  } else if (cens.beg && cens.end) {# Censoring at the beginning and at the end
    
    loglik <- function(par) {
      
      mask <- res$Nik[i, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      fk <- rep.int(x = 0, times = Kmax)
      fk[mask] <- dnbinom(x = kmask, size = par[1], prob = par[2], log = TRUE)
      
      mask <- res$Nebik[i, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = Kmax)
      Fk[mask] <- pnbinom(q = kmask, size = par[1], prob = par[2], lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(res$Nik[i, ] * fk) + sum(res$Nebik[i, ] * Fk)))
    }
    
    CO2 <- constrOptim(
      theta = theta0,
      f = loglik,
      ui = rbind(u0, u1),
      ci = c(c0, c1),
      method = "Nelder-Mead"
    )
    theta <- CO2$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(theta)
  
}