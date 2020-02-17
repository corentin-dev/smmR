.fit.param.fi.nbinom <- function(res, i, Kmax, cens.beg, cens.end) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  xbar <- sum(0:(Kmax - 1) * res$Nik[i, ]) / sum(res$Nik[i, ])
  s2 <- (1 / (sum(res$Nik[i, ]) - 1)) * sum(res$Nik[i, ] * (0:(Kmax - 1) - xbar) ^ 2)
  
  if (xbar >= s2) {
    stop("The negative binomial is not the appropriate distribution for modeling the data")
  }
  
  alphahat <- xbar ^ 2 / (s2 - xbar)
  muhat <- xbar
  
  theta0 <- c(alphahat, muhat)
  
  loglik <- function(par) {
    return(
      -(sum(res$Nik[i, ] * dnbinom(x = 0:(Kmax - 1), size = par[1], mu = par[2], log = TRUE)))
    )
  }
  
  CO2 <- optim(par = theta0, loglik, method = "Nelder-Mead")
  theta0 <- CO2$par
  
  if (!cens.beg && cens.end) {# Censoring at the end
    
    loglik <- function(par) {
      return(
        -(sum(res$Nik[i, ] * dnbinom(x = 0:(Kmax - 1), size = par[1], mu = par[2], log = TRUE)) 
          + sum(res$Neik[i, ] * log(pnbinom(q = 0:(Kmax - 1), size = par[1], mu = par[2], lower.tail = FALSE))))
        )
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Nelder-Mead")
    theta <- CO2$par
    
  } else if (cens.beg && !cens.end) {# Censoring at the beginning
    
    loglik <- function(par) {
      return(
        -(sum(res$Nik[i, ] * dnbinom(x = 0:(Kmax - 1), size = par[1], mu = par[2], log = TRUE)) 
          + sum(res$Nbik[i, ] * log(pnbinom(q = 0:(Kmax - 1), size = par[1], mu = par[2], lower.tail = FALSE))))
        )
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Nelder-Mead")
    theta <- CO2$par
    
  } else if (cens.beg && cens.end) {# Censoring at the beginning and at the end
    
    loglik <- function(par) {
      return(
        -(sum(res$Nik[i, ] * dnbinom(x = 0:(Kmax - 1), size = par[1], mu = par[2], log = TRUE)) 
          + sum(res$Nebik[i, ] * log(pnbinom(q = 0:(Kmax - 1), size = par[1], mu = par[2], lower.tail = FALSE))))
        )
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Nelder-Mead")
    theta <- CO2$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(theta)
  
}