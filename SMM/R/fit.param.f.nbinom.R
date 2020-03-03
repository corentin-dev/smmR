.fit.param.f.nbinom <- function(res, Kmax, cens.beg, cens.end) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  xbar <- sum(0:(Kmax - 1) * res$Nk) / sum(res$Nk)
  s2 <- (1 / (sum(res$Nk) - 1)) * sum(res$Nk * (0:(Kmax - 1) - xbar) ^ 2)
  
  if (xbar >= s2) {
    stop("The negative binomial is not the appropriate distribution for modeling the data")
  }
  
  alphahat <- xbar ^ 2 / (s2 - xbar)
  muhat <- xbar
  
  theta0 <- c(alphahat, muhat)
  
  loglik <- function(par) {
    
    mask <- res$Nk != 0
    kmask <- (0:(Kmax - 1))[mask]
    fk <- rep.int(x = 0, times = Kmax)
    fk[mask] <- dnbinom(x = kmask, size = par[1], mu = par[2], log = TRUE)
    
    return(-(sum(res$Nk * fk)))
  }
  
  # Constraints about the values of the parameters:
  
  # alpha, mu > 0
  u0 <- diag(x = 1, nrow = 2)
  c0 <- c(0, 0)
  
  CO2 <- constrOptim(
    theta = theta0,
    f = loglik,
    ui = u0,
    ci = c0,
    method = "Nelder-Mead"
  )
  theta0 <- CO2$par
  
  if (!cens.beg && cens.end) {# Censoring at the end
    
    loglik <- function(par) {
      
      mask <- res$Nk != 0
      kmask <- (0:(Kmax - 1))[mask]
      fk <- rep.int(x = 0, times = Kmax)
      fk[mask] <- dnbinom(x = kmask, size = par[1], mu = par[2], log = TRUE)
      
      mask <- res$Nek != 0
      kmask <- (0:(Kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = Kmax)
      Fk[mask] <- pnbinom(q = kmask, size = par[1], mu = par[2], lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(res$Nk * fk) + sum(res$Nek * Fk)))
    }
    
    CO2 <- constrOptim(
      theta = theta0,
      f = loglik,
      ui = u0,
      ci = c0,
      method = "Nelder-Mead"
    )
    theta <- CO2$par
    
  } else if (cens.beg && !cens.end) {# Censoring at the beginning
    
    loglik <- function(par) {
      
      mask <- res$Nk != 0
      kmask <- (0:(Kmax - 1))[mask]
      fk <- rep.int(x = 0, times = Kmax)
      fk[mask] <- dnbinom(x = kmask, size = par[1], mu = par[2], log = TRUE)
      
      mask <- res$Nbk != 0
      kmask <- (0:(Kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = Kmax)
      Fk[mask] <- pnbinom(q = kmask, size = par[1], mu = par[2], lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(res$Nk * fk) + sum(res$Nbk * Fk)))
    }
    
    CO2 <- constrOptim(
      theta = theta0,
      f = loglik,
      ui = u0,
      ci = c0,
      method = "Nelder-Mead"
    )
    theta <- CO2$par
    
  } else if (cens.beg && cens.end) {# Censoring at the beginning and at the end
    
    loglik <- function(par) {
      
      mask <- res$Nk != 0
      kmask <- (0:(Kmax - 1))[mask]
      fk <- rep.int(x = 0, times = Kmax)
      fk[mask] <- dnbinom(x = kmask, size = par[1], mu = par[2], log = TRUE)
      
      mask <- res$Nebk != 0
      kmask <- (0:(Kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = Kmax)
      Fk[mask] <- pnbinom(q = kmask, size = par[1], mu = par[2], lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(res$Nk * fk) + sum(res$Nebk * Fk)))
    }
    
    CO2 <- constrOptim(
      theta = theta0,
      f = loglik,
      ui = u0,
      ci = c0,
      method = "Nelder-Mead"
    )
    theta <- CO2$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(theta)
  
}