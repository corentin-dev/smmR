.fit.param.fij.nbinom <- function(res, i, j, Kmax, cens.beg) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  xbar <- sum(0:(Kmax - 1) * res$Nijk[i, j, ]) / sum(res$Nijk[i, j, ])
  s2 <- (1 / (sum(res$Nijk[i, j, ]) - 1)) * sum(res$Nijk[i, j, ] * (0:(Kmax - 1) - xbar) ^ 2)
  
  if (xbar >= s2) {
    stop("The negative binomial is not the appropriate distribution for modeling the data")
  }
  
  alphahat <- xbar ^ 2 / (s2 - xbar)
  muhat <- xbar
  
  theta0 <- c(alphahat, muhat)
  
  loglik <- function(par) {
    
    mask <- res$Nijk[i, j, ] != 0
    kmask <- (0:(Kmax - 1))[mask]
    fk <- rep.int(x = 0, times = Kmax)
    fk[mask] <- dnbinom(x = kmask, size = par[1], mu = par[2], log = TRUE)
    
    return(-(sum(res$Nijk[i, j, ] * fk)))
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
  
  if (cens.beg) {# Censoring at the beginning
    
    loglik <- function(par) {
      
      mask <- res$Nijk[i, j, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      fk <- rep.int(x = 0, times = Kmax)
      fk[mask] <- dnbinom(x = kmask, size = par[1], mu = par[2], log = TRUE)
      
      mask <- res$Nbijk[i, j, ] != 0
      kmask <- (0:(Kmax - 1))[mask]
      Fk <- rep.int(x = 0, times = Kmax)
      Fk[mask] <- pnbinom(q = kmask, size = par[1], mu = par[2], lower.tail = FALSE, log.p = TRUE)
      
      return(-(sum(res$Nijk[i, j, ] * fk) + sum(res$Nbijk[i, j, ] * Fk)))
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