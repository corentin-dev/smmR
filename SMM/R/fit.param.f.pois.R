.fit.param.f.pois <- function(res, Kmax, cens.beg, cens.end) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  theta0 <- sum(0:(Kmax - 1) * res$Nk) / sum(res$Nk)
  
  if (!cens.beg && cens.end) {# Censoring at the end
    
    loglik <- function(par) {
      return(
        -(sum(res$Nk * dpois(x = 0:(Kmax - 1), lambda = par, log = TRUE))
          + sum(res$Nek * log(ppois(q = 0:(Kmax - 1), lambda = par, lower.tail = FALSE))))
      )
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = Kmax - 1)
    theta <- CO2$par
    
  } else if (cens.beg && !cens.end) {# Censoring at the beginning
    
    loglik <- function(par) {
      return(
        -(sum(res$Nk * dpois(x = 0:(Kmax - 1), lambda = par, log = TRUE))
          + sum(res$Nbk * log(ppois(q = 0:(Kmax - 1), lambda = par, lower.tail = FALSE))))
      )
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = Kmax - 1)
    theta <- CO2$par
    
  } else if (cens.beg && cens.end) {# Censoring at the beginning and at the end
    
    loglik <- function(par) {
      return(
        -(sum(res$Nk * dpois(x = 0:(Kmax - 1), lambda = par, log = TRUE))
          + sum(res$Nebk * log(ppois(q = 0:(Kmax - 1), lambda = par, lower.tail = FALSE))))
      )
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = Kmax - 1)
    theta <- CO2$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(c(theta, NA))
  
}