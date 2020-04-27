.fit.param.fj.dweibull <- function(res, j, Kmax, cens.beg) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  theta0 <- suppressWarnings(estdweibull(x = unlist(sapply(1:Kmax, function(x) rep(x, res$Njk[j, x]))), method = "ML", zero = FALSE))

  loglik <- function(par) {
    
    mask <- res$Njk[j, ] != 0
    kmask <- (1:Kmax)[mask]
    fk <- rep.int(x = 0, times = Kmax)
    fk[mask] <- log(ddweibull(x = kmask, q = par[1], beta = par[2], zero = FALSE))
    
    return(-(sum(res$Njk[j, ] * fk)))
  }
  
  # Constraints about the values of the parameters:
  
  # q, beta > 0
  u0 <- diag(x = 1, nrow = 2)
  c0 <- c(0, 0)
  
  # q < 1
  u1 <- matrix(data = c(-1, 0), nrow = 1, ncol = 2)
  c1 <- c(-1)
  
  CO2 <- constrOptim(
    theta = theta0,
    f = loglik,
    ui = rbind(u0, u1),
    ci = c(c0, c1),
    method = "Nelder-Mead"
  )
  theta0 <- CO2$par
  
  if (cens.beg) {# Censoring at the beginning
    
    loglik <- function(par) {
      
      mask <- res$Njk[j, ] != 0
      kmask <- (1:Kmax)[mask]
      fk <- rep.int(x = 0, times = Kmax)
      fk[mask] <- log(ddweibull(x = kmask, q = par[1], beta = par[2], zero = FALSE))
      
      mask <- res$Nbjk[j, ] != 0
      kmask <- (1:Kmax)[mask]
      Fk <- rep.int(x = 0, times = Kmax)
      Fk[mask] <- log(1 - pdweibull(x = kmask, q = par[1], beta = par[2], zero = FALSE))
      
      return(-(sum(res$Njk[j, ] * fk) + sum(res$Nbjk[j, ] * Fk)))
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