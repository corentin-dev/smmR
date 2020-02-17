.fit.param.fi.geom <- function(res, i, Kmax, cens.beg, cens.end) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  theta0 <- sum(res$Nik[i, ]) / sum(1:Kmax * res$Nik[i, ])

  if (!cens.beg && cens.end) {# Censoring at the end
    
    loglik <- function(par) {
      return(
        -(sum(res$Nik[i, ] * dgeom(x = 0:(Kmax - 1), prob = par, log = TRUE))
          + sum(res$Neik[i, ] * log(pgeom(q = 0:(Kmax - 1), prob = par, lower.tail = FALSE))))
        )
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = 1)
    theta <- CO2$par

  } else if (cens.beg && !cens.end) {# censoring at the beginning
    
    loglik <- function(par) {
      return(
        -(sum(res$Nik[i, ] * dgeom(x = 0:(Kmax - 1), prob = par, log = TRUE)) 
          + sum(res$Nbik[i, ] * log(pgeom(q = 0:(Kmax - 1), prob = par, lower.tail = FALSE))))
        )
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = 1)
    theta <- CO2$par
    
  } else if (cens.beg && cens.end) {# censoring at the beginningand at the end
    
    loglik <- function(par) {
      return(
        -(sum(res$Nik[i, ] * dgeom(x = 0:(Kmax - 1), prob = par, log = TRUE)) 
          + sum(res$Nebik[i, ] * log(pgeom(q = 0:(Kmax - 1), prob = par, lower.tail = FALSE))))
        )
    }

    CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = 1)
    theta <- CO2$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(c(theta, NA))
    
}