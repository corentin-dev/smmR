.fit.param.f.dweibull <- function(res, Kmax, cens.beg, cens.end) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  theta0 <- suppressMessages(estdweibull(x = unlist(sapply(1:Kmax, function(x) rep(x, res$Nk[x]))), method = "ML", zero = FALSE))
  
  loglik <- function(par) {
    return(
      -(sum(res$Nk * log(ddweibull(x = 1:Kmax, q = par[1], beta = par[2], zero = FALSE))))
    )
  }

  CO2 <- optim(par = theta0, loglik, method = "Nelder-Mead")
  theta0 <- CO2$par
  
  if (!cens.beg && cens.end) {# Censoring at the end
    
    loglik <- function(par) {
      return(
        -(sum(res$Nk * log(ddweibull(x = 1:Kmax, q = par[1], beta = par[2], zero = FALSE))) 
          + sum(res$Nek * log(1 - pdweibull(x = 1:Kmax, q = par[1], beta = par[2], zero = FALSE))))
        )
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Nelder-Mead")
    theta <- CO2$par
    
  } else if (cens.beg && !cens.end) {# Censoring at the beginning
    
    loglik <- function(par) {
      return(
        -(sum(res$Nk * log(ddweibull(x = 1:Kmax, q = par[1], beta = par[2], zero = FALSE))) 
          + sum(res$Nbk * log(1 - pdweibull(x = 1:Kmax, q = par[1], beta = par[2], zero = FALSE))))
        )
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Nelder-Mead")
    theta <- CO2$par
    
  } else if (cens.beg && cens.end) {# Censoring at the beginning and at the end
    
    loglik <- function(par) {
      return(
        -(sum(res$Nk * log(ddweibull(x = 1:Kmax, q = par[1], beta = par[2], zero = FALSE))) 
          + sum(res$Nebk * log(1 - pdweibull(x = 1:Kmax, q = par[1], beta = par[2], zero = FALSE))))
        )
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Nelder-Mead")
    theta <- CO2$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(theta)
  
}