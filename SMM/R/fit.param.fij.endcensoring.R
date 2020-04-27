.fit.param.fij.endcensoring <- function(res, i, S, Kmax, distr, cens.beg) {
  
  # We use the non-censoring case as initial values for the optimization problem
  theta0 <- matrix(data = NA, nrow = 2, ncol = S)
  for (j in 1:S) {
    if (j != i) {
      if (distr[i, j] == "dweibull") {
        theta0[, j] <- .fit.param.fij.dweibull(res, i, j, Kmax, cens.beg = FALSE)
      } else if (distr[i, j] == "geom") {
        theta0[, j] <- .fit.param.fij.geom(res, i, j, Kmax, cens.beg = FALSE)
      } else if (distr[i, j] == "nbinom") {
        theta0[, j] <- .fit.param.fij.nbinom(res, i, j, Kmax, cens.beg = FALSE)
      } else if (distr[i, j] == "pois") {
        theta0[, j] <- .fit.param.fij.pois(res, i, j, Kmax, cens.beg = FALSE)
      }
    }
  }
  
  theta0 <- as.vector(theta0[!is.na(theta0)])
  
  ##########################################################
  # Specific case S = 2
  ##########################################################
  if (S == 2) {
    
    if (i == 1) {
      puv <- c(0, 1)
    } else if (i == 2) {
      puv <- c(1, 0)
    }
    
    param <- matrix(data = NA, nrow = S, ncol = 2)
    
    ##########################################################
    # Sojourn time distribution is a uniform
    ##########################################################
    if ("unif" %in% distr[i, ]) {
      
      param[which(distr[i, ] == "unif"), ] <- .fit.param.fij.unif(res, i, which(distr[i, ] == "unif"), Kmax)
      
    } else {
      
      if (cens.beg) {
        
        loglik <- function(par) {
          
          fv <- matrix(data = 0, nrow = S, ncol = Kmax)
          Fv <- matrix(data = 0, nrow = S, ncol = Kmax)
          fv2 <- matrix(data = 0, nrow = S, ncol = Kmax)
          
          maskNeik <- res$Neik[i, ] != 0
          
          skipindex <- 1
          for (j in 1:S) {
            
            if (is.na(distr[i, j])) {
              
              fv[j, ] <- rep.int(x = 0, times = Kmax)
              Fv[j, ] <- rep.int(x = 0, times = Kmax)
              fv2[j, ] <- rep.int(x = 0, times = Kmax)
              
            } else {
              
              maskNijk <- res$Nijk[i, j, ] != 0
              maskNbijk <- res$Nbijk[i, j, ] != 0
              
              if (distr[i, j] == "dweibull") {
                fv[j, maskNijk] <- log(ddweibull(x = (1:Kmax)[maskNijk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
                Fv[j, maskNbijk] <- log(1 - pdweibull(x = (1:Kmax)[maskNbijk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
                fv2[j, ] <- ddweibull(x = 1:Kmax, q = par[skipindex], beta = par[skipindex + 1], zero = FALSE)
                skipindex <- skipindex + 2
              } else if (distr[i, j] == "geom") {
                fv[j, maskNijk] <- dgeom(x = (0:(Kmax - 1))[maskNijk], prob = par[skipindex], log = TRUE)
                Fv[j, maskNbijk] <- pgeom(q = (0:(Kmax - 1))[maskNbijk], prob = par[skipindex], lower.tail = FALSE, log.p = TRUE)
                fv2[j, ] <- dgeom(x = 0:(Kmax - 1), prob = par[skipindex])
                skipindex <- skipindex + 1
              } else if (distr[i, j] == "nbinom") {
                fv[j, maskNijk] <- dnbinom(x = (0:(Kmax - 1))[maskNijk], size = par[skipindex], prob = par[skipindex + 1], log = TRUE)
                Fv[j, maskNbijk] <- pnbinom(q = (0:(Kmax - 1))[maskNbijk], size = par[skipindex], prob = par[skipindex + 1], lower.tail = FALSE, log.p = TRUE)
                fv2[j, ] <- dnbinom(x = 0:(Kmax - 1), size = par[skipindex], prob = par[skipindex + 1])
                skipindex <- skipindex + 2
              } else if (distr[i, j] == "pois") {
                fv[j, maskNijk] <- dpois(x = (0:(Kmax - 1))[maskNijk], lambda = par[skipindex], log = TRUE)
                Fv[j, maskNbijk] <- ppois(q = (0:(Kmax - 1))[maskNbijk], lambda = par[skipindex], lower.tail = FALSE, log.p = TRUE)
                fv2[j, ] <- dpois(x = 0:(Kmax - 1), lambda = par[skipindex])
                skipindex <- skipindex + 1
              }
            }
          }
          
          fv2 <- apply(apply(puv * fv2, 1, cumsum), 1, sum)
          fv2[maskNeik] <- log(1 - fv2[maskNeik])
          
          return(
            -(sum(res$Nijk[i, , ] * fv) + sum(res$Nbijk[i, , ] * Fv)
              + sum(res$Neik[i, ] * fv2))
          )
        }
        
      } else {
        
        loglik <- function(par) {
          
          fv <- matrix(data = 0, nrow = S, ncol = Kmax)
          fv2 <- matrix(data = 0, nrow = S, ncol = Kmax)
          
          maskNeik <- res$Neik[i, ] != 0
          
          skipindex <- 1
          for (j in 1:S) {
            if (is.na(distr[i, j])) {
              
              fv[j, ] <- rep.int(x = 0, times = Kmax)
              fv2[j, ] <- rep.int(x = 0, times = Kmax)
              
            } else {
              
              maskNijk <- res$Nijk[i, j, ] != 0
              
              if (distr[i, j] == "dweibull") {
                fv[j, maskNijk] <- log(ddweibull(x = (1:Kmax)[maskNijk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
                fv2[j, ] <- ddweibull(x = 1:Kmax, q = par[skipindex], beta = par[skipindex + 1], zero = FALSE)
                skipindex <- skipindex + 2
              } else if (distr[i, j] == "geom") {
                fv[j, maskNijk] <- dgeom(x = (0:(Kmax - 1))[maskNijk], prob = par[skipindex], log = TRUE)
                fv2[j, ] <- dgeom(x = 0:(Kmax - 1), prob = par[skipindex])
                skipindex <- skipindex + 1
              } else if (distr[i, j] == "nbinom") {
                fv[j, maskNijk] <- dnbinom(x = (0:(Kmax - 1))[maskNijk], size = par[skipindex], prob = par[skipindex + 1], log = TRUE)
                fv2[j, ] <- dnbinom(x = 0:(Kmax - 1), size = par[skipindex], prob = par[skipindex + 1])
                skipindex <- skipindex + 2
              } else if (distr[i, j] == "pois") {
                fv[j, maskNijk] <- dpois(x = (0:(Kmax - 1))[maskNijk], lambda = par[skipindex], log = TRUE)
                fv2[j, ] <- dpois(x = 0:(Kmax - 1), lambda = par[skipindex])
                skipindex <- skipindex + 1
              }
            }
          }
          
          fv2 <- apply(apply(puv * fv2, 1, cumsum), 1, sum)
          fv2[maskNeik] <- log(1 - fv2[maskNeik])
          
          return(
            -(sum(res$Nijk[i, , ] * fv) + sum(res$Neik[i, ] * fv2))
          )
        }
        
      }
      
      if (distr[i, abs(i - 3)] == "dweibull") {
        
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
        param[abs(i - 3), ] <- CO2$par
        
      } else if (distr[i, abs(i - 3)] == "geom") {
        
        CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = 1)
        param[abs(i - 3), ] <- CO2$par
        
      } else if (distr[i, abs(i - 3)] == "nbinom") {
        
        # Constraints about the values of the parameters
        
        # alpha, p > 0
        u0 <- diag(x = 1, nrow = 2)
        c0 <- c(0, 0)
        
        # p < 1
        u1 <- matrix(data = c(0, -1), nrow = 1, ncol = 2)
        c1 <- c(-1)
        
        CO2 <- constrOptim(
          theta = theta0,
          f = loglik,
          ui = rbind(u0, u1),
          ci = c(c0, c1),
          method = "Nelder-Mead"
        )
        param[abs(i - 3), ] <- CO2$par
        
      } else if (distr[i, abs(i - 3)] == "pois") {
        
        CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = Kmax - 1)
        param[abs(i - 3), ] <- CO2$par
        
      }
      
    }
    
    ##########################################################
    # Case S > 2
    ########################################################## 
  } else {
    
    if (cens.beg) {
      
      loglik <- function(par) {
        
        parpuv <- par[1:(S - 2)]
        parpuv <- c(parpuv, 1 - sum(parpuv))
        
        fv <- matrix(data = 0, nrow = S, ncol = Kmax)
        Fv <- matrix(data = 0, nrow = S, ncol = Kmax)
        fv2 <- matrix(data = 0, nrow = S, ncol = Kmax)
        
        maskNeik <- res$Neik[i, ] != 0
        
        skipindex <- (S - 2) + 1
        for (j in 1:S) {
          
          if (is.na(distr[i, j])) {
            
            fv[j, ] <- rep.int(x = 0, times = Kmax)
            Fv[j, ] <- rep.int(x = 0, times = Kmax)
            fv2[j, ] <- rep.int(x = 0, times = Kmax)
            
          } else {
            
            maskNijk <- res$Nijk[i, j, ] != 0
            maskNbijk <- res$Nbijk[i, j, ] != 0
            
            if (distr[i, j] == "dweibull") {
              fv[j, maskNijk] <- log(ddweibull(x = (1:Kmax)[maskNijk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
              Fv[j, maskNbijk] <- log(1 - pdweibull(x = (1:Kmax)[maskNbijk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
              fv2[j, ] <- ddweibull(x = 1:Kmax, q = par[skipindex], beta = par[skipindex + 1], zero = FALSE)
              skipindex <- skipindex + 2
            } else if (distr[i, j] == "geom") {
              fv[j, maskNijk] <- dgeom(x = (0:(Kmax - 1))[maskNijk], prob = par[skipindex], log = TRUE)
              Fv[j, maskNbijk] <- log(pgeom(q = (0:(Kmax - 1))[maskNbijk], prob = par[skipindex], lower.tail = FALSE))
              fv2[j, ] <- dgeom(x = 0:(Kmax - 1), prob = par[skipindex])
              skipindex <- skipindex + 1
            } else if (distr[i, j] == "nbinom") {
              fv[j, maskNijk] <- dnbinom(x = (0:(Kmax - 1))[maskNijk], size = par[skipindex], prob = par[skipindex + 1], log = TRUE)
              Fv[j, maskNbijk] <- pnbinom(q = (0:(Kmax - 1))[maskNbijk], size = par[skipindex], prob = par[skipindex + 1], lower.tail = FALSE, log.p = TRUE)
              fv2[j, ] <- dnbinom(x = 0:(Kmax - 1), size = par[skipindex], prob = par[skipindex + 1])
              skipindex <- skipindex + 2
            } else if (distr[i, j] == "pois") {
              fv[j, maskNijk] <- dpois(x = (0:(Kmax - 1))[maskNijk], lambda = par[skipindex], log = TRUE)
              Fv[j, maskNbijk] <- ppois(q = (0:(Kmax - 1))[maskNbijk], lambda = par[skipindex], lower.tail = FALSE, log.p = TRUE)
              fv2[j, ] <- dpois(x = 0:(Kmax - 1), lambda = par[skipindex])
              skipindex <- skipindex + 1
            }
          }
        }
        
        fv2 <- apply(apply(parpuv * fv2[-i, ], 1, cumsum), 1, sum)
        fv2[maskNeik] <- log(1 - fv2[maskNeik])
        
        return(-(
          sum(res$Nij[i, -i] * log(parpuv)) +
            sum(res$Nijk[i, , ] * fv) + sum(res$Nbijk[i, , ] * Fv) + 
            sum(res$Neik[i, ] * fv2)
        ))
      }
      
    } else {
      
      loglik <- function(par) {
        
        parpuv <- par[1:(S - 2)]
        parpuv <- c(parpuv, 1 - sum(parpuv))
        
        fv <- matrix(data = 0, nrow = S, ncol = Kmax)
        fv2 <- matrix(data = 0, nrow = S, ncol = Kmax)
        
        maskNeik <- res$Neik[i, ] != 0
        
        skipindex <- (S - 2) + 1
        for (j in 1:S) {
          
          if (is.na(distr[i, j])) {
            
            fv[j, ] <- rep.int(x = 0, times = Kmax)
            fv2[j, ] <- rep.int(x = 0, times = Kmax)
            
          } else {
            
            maskNijk <- res$Nijk[i, j, ] != 0
            
            if (distr[i, j] == "dweibull") {
              fv[j, maskNijk] <- log(ddweibull(x = (1:Kmax)[maskNijk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
              fv2[j, ] <- ddweibull(x = 1:Kmax, q = par[skipindex], beta = par[skipindex + 1], zero = FALSE)
              skipindex <- skipindex + 2
            } else if (distr[i, j] == "geom") {
              fv[j, maskNijk] <- dgeom(x = (0:(Kmax - 1))[maskNijk], prob = par[skipindex], log = TRUE)
              fv2[j, ] <- dgeom(x = 0:(Kmax - 1), prob = par[skipindex])
              skipindex <- skipindex + 1
            } else if (distr[i, j] == "nbinom") {
              fv[j, maskNijk] <- dnbinom(x = (0:(Kmax - 1))[maskNijk], size = par[skipindex], prob = par[skipindex + 1], log = TRUE)
              fv2[j, ] <- dnbinom(x = 0:(Kmax - 1), size = par[skipindex], prob = par[skipindex + 1])
              skipindex <- skipindex + 2
            } else if (distr[i, j] == "pois") {
              fv[j, maskNijk] <- dpois(x = (0:(Kmax - 1))[maskNijk], lambda = par[skipindex], log = TRUE)
              fv2[j, ] <- dpois(x = 0:(Kmax - 1), lambda = par[skipindex])
              skipindex <- skipindex + 1
            }
          }
        }
        
        fv2 <- apply(apply(parpuv * fv2[-i, ], 1, cumsum), 1, sum)
        fv2[maskNeik] <- log(1 - fv2[maskNeik])
        
        return(-(
          sum(res$Nij[i, -i] * log(parpuv)) +
            sum(res$Nijk[i, , ] * fv) + 
            sum(res$Neik[i, ] * fv2)
        ))
      }
      
    }
    
    # Constraints about the values of the parameters:
    
    # puv >= 0 and parameters values > 0
    u0 <- diag(x = 1, nrow = (S - 2) + length(theta0), ncol = (S - 2) + length(theta0))
    c0 <- rep(0, (S - 2) + length(theta0))
    
    # puv <= 1
    u1 <- matrix(0, nrow = (S - 2), ncol = (S - 2) + length(theta0))
    diag(u1) <- -1
    c1 <- rep(-1, (S - 2))
    
    # sum(puv) <= 1
    u2 <- matrix(0, nrow = 1, ncol = (S - 2) + length(theta0))
    u2[1, 1:(S - 2)] <- -1
    c2 <- c(-1)
    
    # Specific constraints depending on the distribution
    u3 <- matrix(0, nrow = length(theta0), ncol = (S - 2) + length(theta0))
    skipindex <- 1
    rowstoremove <- c()
    for (j in 1:S) {
      if (j != i) {
        if (distr[i, j] == "dweibull") {
          u3[skipindex, skipindex + (S - 2)] <- -1
          rowstoremove <- c(rowstoremove, skipindex + 1)
          skipindex <- skipindex + 2
        } else if (distr[i, j] == "geom") {
          u3[skipindex, skipindex + (S - 2)] <- -1
          skipindex <- skipindex + 1
        } else if (distr[i, j] == "nbinom") {
          u3[skipindex + 1, skipindex + 1 + (S - 2)] <- -1
          rowstoremove <- c(rowstoremove, skipindex)
          skipindex <- skipindex + 2
        } else if (distr[i, j] == "pois") {
          rowstoremove <- c(rowstoremove, skipindex)
          skipindex <- skipindex + 1
        }
      }
    }
    
    if (!is.null(rowstoremove)) {
      u3 <- u3[-rowstoremove, ] 
    }
    
    if (!is.null(nrow(u3))) {
      c3 <- rep(-1, nrow(u3))
    } else {
      c3 <- c(-1)
    }
    
    if (length(u3) != 0) {
      u4 <- rbind(u0, u1, u2, u3)
      c4 <- c(c0, c1, c2, c3)
    } else {
      u4 <- rbind(u0, u1, u2)
      c4 <- c(c0, c1, c2)
    }
    
    CO2 <-
      constrOptim(
        theta = c(rep(1 / (S - 1), (S - 2)), theta0),
        loglik,
        ui = u4,
        ci = c4,
        method = "Nelder-Mead"
      )
    
    # Rebuild parameters
    parpuv <- CO2$par[1:(S - 2)]
    parpuv <- c(parpuv, 1 - sum(parpuv))
    
    if (i == 1) {
      parpuv <- c(0, parpuv)
    } else if (i == S) {
      parpuv <- c(parpuv, 0)
    } else {
      parpuv <- c(parpuv[1:(i - 1)], 0, parpuv[i:length(parpuv)])
    }
    
    param <- matrix(data = NA, nrow = S, ncol = 2)
    skipindex <- (S - 2) + 1
    for (j in 1:S) {
      if (j != i) {
        if (distr[i, j] %in% c("dweibull", "nbinom")) {
          param[j, ] <- CO2$par[skipindex:(skipindex + 1)]
          skipindex <- skipindex + 2
        } else if (distr[i, j] %in% c("geom", "pois")) {
          param[j, 1] <- CO2$par[skipindex]
          skipindex <- skipindex + 1
        } else if (distr[i, j] == "unif") {
          param[j, ] <- .fit.param.fij.unif(res, i, j, Kmax)
        }
      }
    }
    
    puv <- parpuv
    
  }
  
  return(list(ptrans = puv, param = param))
}