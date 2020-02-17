.fit.param.fi.dweibull <- function(res, i, Kmax, cens.beg, cens.end) {
  
  # Estimation of the parameters of the distribution (No censoring case)
  # try(theta0 <- suppressWarnings(suppressMessages(estdweibull(x = unlist(sapply(1:Kmax, function(x) rep(x, res$Nik[i, x]))), method = "M", zero = FALSE))), silent = TRUE)
  try(theta0 <- suppressWarnings(suppressMessages(estdweibull(x = unlist(sapply(1:Kmax, function(x) rep(x, res$Nik[i, x]))), method = "P", zero = FALSE))), silent = TRUE)
  try(theta0 <- suppressWarnings(suppressMessages(estdweibull(x = unlist(sapply(1:Kmax, function(x) rep(x, res$Nik[i, x]))), method = "ML", zero = FALSE))), silent = TRUE)

  loglik <- function(par) {
    return(
      -(sum(res$Nik[i, ] * log(ddweibull(x = 1:Kmax, q = par[1], beta = par[2], zero = FALSE))))
    )
  }
  
  CO2 <- optim(par = theta0, loglik, method = "Nelder-Mead")
  theta0 <- CO2$par
  
  if (!cens.beg && cens.end) {# Censoring at the end
    
    loglik <- function(par) {
      return(
        -(sum(res$Nik[i, ] * log(ddweibull(x = 1:Kmax, q = par[1], beta = par[2], zero = FALSE))) 
          + sum(res$Neik[i, ] * log(1 - pdweibull(x = 1:Kmax, q = par[1], beta = par[2], zero = FALSE))))
        )
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Nelder-Mead")
    theta <- CO2$par
    
  } else if (cens.beg && !cens.end) {# Censoring at the beginning
    
    loglik <- function(par) {
      return(
        -(sum(res$Nik[i, ] * log(ddweibull(x = 1:Kmax, q = par[1], beta = par[2], zero = FALSE))) 
          + sum(res$Nbik[i, ] * log(1 - pdweibull(x = 1:Kmax, q = par[1], beta = par[2], zero = FALSE))))
        )
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Nelder-Mead")
    theta <- CO2$par
    
  } else if (cens.beg && cens.end) {# Censoring at the beginning and at the end
    
    loglik <- function(par) {
      return(
        -(sum(res$Nik[i, ] * log(ddweibull(x = 1:Kmax, q = par[1], beta = par[2], zero = FALSE))) 
          + sum(res$Nebik[i, ] * log(1 - pdweibull(x = 1:Kmax, q = par[1], beta = par[2], zero = FALSE))))
        )
    }
    
    CO2 <- optim(par = theta0, loglik, method = "Nelder-Mead")
    theta <- CO2$par
    
  } else {# No censoring
    
    theta <- theta0
    
  }
  
  return(theta)
  
}