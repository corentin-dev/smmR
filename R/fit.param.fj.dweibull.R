.fit.param.fj.dweibull <- function(counting, j, kmax, cens.beg) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  theta0 <- suppressWarnings(estdweibull(x = unlist(sapply(1:kmax, function(x) rep(x, counting$Njk[j, x]))), method = "ML", zero = FALSE))
  
  logLik <- function(par) {
    
    mask <- counting$Njk[j, ] != 0
    kmask <- (1:kmax)[mask]
    fk <- rep.int(x = 0, times = kmax)
    fk[mask] <- log(ddweibull(x = kmask, q = par[1], beta = par[2], zero = FALSE))
    
    return(-(sum(counting$Njk[j, ] * fk)))
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
    f = logLik,
    ui = rbind(u0, u1),
    ci = c(c0, c1),
    method = "Nelder-Mead"
  )
  theta0 <- mle$par
  
  if (cens.beg) {# Censoring at the beginning
    
    logLik <- function(par) {
      
      mask <- counting$Njk[j, ] != 0
      kmask <- (1:kmax)[mask]
      fk <- rep.int(x = 0, times = kmax)
      fk[mask] <- log(ddweibull(x = kmask, q = par[1], beta = par[2], zero = FALSE))
      
      mask <- counting$Nbjk[j, ] != 0
      kmask <- (1:kmax)[mask]
      Fk <- rep.int(x = 0, times = kmax)
      Fk[mask] <- log(1 - pdweibull(x = kmask, q = par[1], beta = par[2], zero = FALSE))
      
      return(-(sum(counting$Njk[j, ] * fk) + sum(counting$Nbjk[j, ] * Fk)))
    }
    
    mle <- constrOptim(
      theta = theta0,
      f = logLik,
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
