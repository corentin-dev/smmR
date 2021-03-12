.fit.param.f.nbinom <- function(counting, kmax, cens.beg, cens.end) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  xbar <- sum(0:(kmax - 1) * counting$Nk) / sum(counting$Nk)
  s2 <- (1 / (sum(counting$Nk) - 1)) * sum(counting$Nk * (0:(kmax - 1) - xbar) ^ 2)
  
  if (xbar >= s2) {
    stop("The negative binomial distribution is not appropriate for 
         modeling the sojourn time distribution (variance < expectation)")
  }
  
  alphahat <- xbar ^ 2 / (s2 - xbar)
  phat <- xbar / s2
  
  theta0 <- c(alphahat, phat)
  
  loglik <- function(par) {
    
    mask <- counting$Nk != 0
    kmask <- (0:(kmax - 1))[mask]
    fk <- rep.int(x = 0, times = kmax)
    fk[mask] <- dnbinom(x = kmask, size = par[1], prob = par[2], log = TRUE)
    
    return(-(sum(counting$Nk * fk)))
  }
  
  # Constraints about the values of the parameters:
  
  # alpha, p > 0
  u0 <- diag(x = 1, nrow = 2)
  c0 <- c(0, 0)
  
  # p < 1
  u1 <- matrix(data = c(0, -1), nrow = 1, ncol = 2)
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
      
      mask <- counting$Nk != 0
      kmask <- (0:(kmax - 1))[mask]
      fk <- rep.int(x = 0, times = kmax)
      fk[mask] <- dnbinom(x = kmask, size = par[1], prob = par[2], log = TRUE)
      
      mask <- counting$Nek != 0
      kmask <- (0:(kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = kmax)
      Fk[mask] <- pnbinom(q = kmask, size = par[1], prob = par[2], lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(counting$Nk * fk) + sum(counting$Nek * Fk)))
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
      
      mask <- counting$Nk != 0
      kmask <- (0:(kmax - 1))[mask]
      fk <- rep.int(x = 0, times = kmax)
      fk[mask] <- dnbinom(x = kmask, size = par[1], prob = par[2], log = TRUE)
      
      mask <- counting$Nbk != 0
      kmask <- (0:(kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = kmax)
      Fk[mask] <- pnbinom(q = kmask, size = par[1], prob = par[2], lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(counting$Nk * fk) + sum(counting$Nbk * Fk)))
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
      
      mask <- counting$Nk != 0
      kmask <- (0:(kmax - 1))[mask]
      fk <- rep.int(x = 0, times = kmax)
      fk[mask] <- dnbinom(x = kmask, size = par[1], prob = par[2], log = TRUE)
      
      mask <- counting$Nebk != 0
      kmask <- (0:(kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = kmax)
      Fk[mask] <- pnbinom(q = kmask, size = par[1], prob = par[2], lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(counting$Nk * fk) + sum(counting$Nebk * Fk)))
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
