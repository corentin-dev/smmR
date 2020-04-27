.fit.param.fj.endcensoring <- function(res, S, Kmax, distr, cens.beg) {
  
  # We use the non-censoring case as initial values for the optimization problem
  theta0 <- matrix(data = NA, nrow = 2, ncol = S)
  for (j in 1:S) {
    
    if (distr[j] == "dweibull") {
      theta0[, j] <- .fit.param.fj.dweibull(res, j, Kmax, cens.beg = FALSE)
    } else if (distr[j] == "geom") {
      theta0[, j] <- .fit.param.fj.geom(res, j, Kmax, cens.beg = FALSE)
    } else if (distr[j] == "nbinom") {
      theta0[, j] <- .fit.param.fj.nbinom(res, j, Kmax, cens.beg = FALSE)
    } else if (distr[j] == "pois") {
      theta0[, j] <- .fit.param.fj.pois(res, j, Kmax, cens.beg = FALSE)
    }
    
  }
  
  theta0 <- as.vector(theta0[!(is.na(theta0))])
  
  ##########################################################
  # Specific case S = 2
  ##########################################################
  if (S == 2) {
    
    ptrans <- matrix(c(0, 1, 1, 0), nrow = S, byrow = TRUE)
    param <- matrix(data = NA, nrow = S, ncol = 2)
    
    ##########################################################
    # Sojourn time distributions are 2 uniforms
    ##########################################################
    if ((distr[1] == "unif") && (distr[2] == "unif")) {
      
      for (j in 1:S) {
        param[j, ] <- .fit.param.fj.unif(res, j)
      }
      
      ##########################################################
      # There is one uniform distribution among the 2 distributions
      ##########################################################
    } else if ("unif" %in% distr) {
      
      for (j in 1:S) {
        if (distr[j] == "unif") {
          param[j, ] <- .fit.param.fj.unif(res, j)  
        } else {
          
          if (cens.beg) {
            
            loglik <- function(par) {
              
              fv <- matrix(data = 0, nrow = S, ncol = Kmax)
              Fv <- matrix(data = 0, nrow = S, ncol = Kmax)
              skipindex <- 1
              for (j in 1:S) {
                
                maskNjk <- res$Njk[j, ] != 0
                maskNeNb <- (res$Nbjk[j, ] + res$Neik[abs(j - 3), ]) != 0
                
                if (distr[j] == "dweibull") {
                  fv[j, maskNjk] <- log(ddweibull(x = (1:Kmax)[maskNjk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
                  Fv[j, maskNeNb] <- log(1 - pdweibull(x = (1:Kmax)[maskNeNb], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE) + .Machine$double.xmin)
                  skipindex <- skipindex + 2
                } else if (distr[j] == "geom") {
                  fv[j, maskNjk] <- dgeom(x = (0:(Kmax - 1))[maskNjk], prob = par[skipindex], log = TRUE)
                  Fv[j, maskNeNb] <- pgeom(q = (0:(Kmax - 1))[maskNeNb], prob = par[skipindex], lower.tail = FALSE, log.p = TRUE)
                  skipindex <- skipindex + 1
                } else if (distr[j] == "nbinom") {
                  fv[j, maskNjk] <- dnbinom(x = (0:(Kmax - 1))[maskNjk], size = par[skipindex], prob = par[skipindex + 1], log = TRUE)
                  Fv[j, maskNeNb] <- pnbinom(q = (0:(Kmax - 1))[maskNeNb], size = par[skipindex], prob = par[skipindex + 1], lower.tail = FALSE, log.p = TRUE)
                  skipindex <- skipindex + 2
                } else if (distr[j] == "pois") {
                  fv[j, maskNjk] <- dpois(x = (0:(Kmax - 1))[maskNjk], lambda = par[skipindex], log = TRUE)
                  Fv[j, maskNeNb] <- ppois(q = (0:(Kmax - 1))[maskNeNb], lambda = par[skipindex], lower.tail = FALSE, log.p = TRUE)
                  skipindex <- skipindex + 1
                }
              }
              
              return(
                -(sum(res$Njk[1, ] * fv[1, ]) + sum(res$Njk[2, ] * fv[2, ])
                  + sum((res$Nbjk[1, ] + res$Neik[2, ]) * Fv[1, ]) 
                  + sum((res$Nbjk[2, ] + res$Neik[1, ]) * Fv[2, ]))
              )
            }
            
          } else {
            
            loglik <- function(par) {
              
              fv <- matrix(data = 0, nrow = S, ncol = Kmax)
              Fv <- matrix(data = 0, nrow = S, ncol = Kmax)
              skipindex <- 1
              for (j in 1:S) {
                
                maskNjk <- res$Njk[j, ] != 0
                maskNeik <- res$Neik[abs(j - 3), ] != 0
                
                if (distr[j] == "dweibull") {
                  fv[j, maskNjk] <- log(ddweibull(x = (1:Kmax)[maskNjk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
                  Fv[j, maskNeik] <- log(1 - pdweibull(x = (1:Kmax)[maskNeik], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE) + .Machine$double.xmin)
                  skipindex <- skipindex + 2
                } else if (distr[j] == "geom") {
                  fv[j, maskNjk] <- dgeom(x = (0:(Kmax - 1))[maskNjk], prob = par[skipindex], log = TRUE)
                  Fv[j, maskNeik] <- pgeom(q = (0:(Kmax - 1))[maskNeik], prob = par[skipindex], lower.tail = FALSE, log.p = TRUE)
                  skipindex <- skipindex + 1
                } else if (distr[j] == "nbinom") {
                  fv[j, maskNjk] <- dnbinom(x = (0:(Kmax - 1))[maskNjk], size = par[skipindex], prob = par[skipindex + 1], log = TRUE)
                  Fv[j, maskNeik] <- pnbinom(q = (0:(Kmax - 1))[maskNeik], size = par[skipindex], prob = par[skipindex + 1], lower.tail = FALSE, log.p = TRUE)
                  skipindex <- skipindex + 2
                } else if (distr[j] == "pois") {
                  fv[j, maskNjk] <- dpois(x = (0:(Kmax - 1))[maskNjk], lambda = par[skipindex], log = TRUE)
                  Fv[j, maskNeik] <- ppois(q = (0:(Kmax - 1))[maskNeik], lambda = par[skipindex], lower.tail = FALSE, log.p = TRUE)
                  skipindex <- skipindex + 1
                }
              }
              
              return(
                -(sum(res$Njk[1, ] * fv[1, ]) + sum(res$Njk[2, ] * fv[2, ])
                  + sum(res$Neik[1, ] * Fv[2, ]) + sum(res$Neik[2, ] * Fv[1, ]))
              )
            }
            
          }
          
          if (distr[j] == "dweibull") {
            
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
            
            param[j, ] <- CO2$par
            
          } else if (distr[j] == "geom") {
            
            CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = 1)
            param[j, ] <- CO2$par
            
          } else if (distr[j] == "nbinom") {
            
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
            param[j, ] <- CO2$par
            
          } else if (distr[j] == "pois") {
            
            CO2 <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = Kmax - 1)
            param[j, ] <- CO2$par
            
          }
        }
      }
      
      ##########################################################
      # Two sojourn time distributions different from the uniform
      ##########################################################
    } else {
      
      if (cens.beg) {
        
        loglik <- function(par) {
          
          fv <- matrix(data = 0, nrow = S, ncol = Kmax)
          Fv <- matrix(data = 0, nrow = S, ncol = Kmax)
          skipindex <- 1
          for (j in 1:S) {
            
            maskNjk <- res$Njk[j, ] != 0
            maskNeNb <- (res$Nbjk[j, ] + res$Neik[abs(j - 3), ]) != 0
            
            if (distr[j] == "dweibull") {
              fv[j, maskNjk] <- log(ddweibull(x = (1:Kmax)[maskNjk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
              Fv[j, maskNeNb] <- log(1 - pdweibull(x = (1:Kmax)[maskNeNb], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE) + .Machine$double.xmin)
              skipindex <- skipindex + 2
            } else if (distr[j] == "geom") {
              fv[j, maskNjk] <- dgeom(x = (0:(Kmax - 1))[maskNjk], prob = par[skipindex], log = TRUE)
              Fv[j, maskNeNb] <- pgeom(q = (0:(Kmax - 1))[maskNeNb], prob = par[skipindex], lower.tail = FALSE, log.p = TRUE)
              skipindex <- skipindex + 1
            } else if (distr[j] == "nbinom") {
              fv[j, maskNjk] <- dnbinom(x = (0:(Kmax - 1))[maskNjk], size = par[skipindex], prob = par[skipindex + 1], log = TRUE)
              Fv[j, maskNeNb] <- pnbinom(q = (0:(Kmax - 1))[maskNeNb], size = par[skipindex], prob = par[skipindex + 1], lower.tail = FALSE, log.p = TRUE)
              skipindex <- skipindex + 2
            } else if (distr[j] == "pois") {
              fv[j, maskNjk] <- dpois(x = (0:(Kmax - 1))[maskNjk], lambda = par[skipindex], log = TRUE)
              Fv[j, maskNeNb] <- ppois(q = (0:(Kmax - 1))[maskNeNb], lambda = par[skipindex], lower.tail = FALSE, log.p = TRUE)
              skipindex <- skipindex + 1
            }
          }
          
          return(
            -(sum(res$Njk[1, ] * fv[1, ]) + sum(res$Njk[2, ] * fv[2, ])
              + sum((res$Nbjk[1, ] + res$Neik[2, ]) * Fv[1, ]) 
              + sum((res$Nbjk[2, ] + res$Neik[1, ]) * Fv[2, ])
            )
          )
        }
        
      } else {
        
        loglik <- function(par) {
          
          fv <- matrix(data = 0, nrow = S, ncol = Kmax)
          Fv <- matrix(data = 0, nrow = S, ncol = Kmax)
          skipindex <- 1
          for (j in 1:S) {
            
            maskNjk <- res$Njk[j, ] != 0
            maskNeik <- res$Neik[abs(j - 3), ] != 0
            
            if (distr[j] == "dweibull") {
              fv[j, maskNjk] <- log(ddweibull(x = (1:Kmax)[maskNjk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
              Fv[j, maskNeik] <- log(1 - pdweibull(x = (1:Kmax)[maskNeik], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE) + .Machine$double.xmin)
              skipindex <- skipindex + 2
            } else if (distr[j] == "geom") {
              fv[j, maskNjk] <- dgeom(x = (0:(Kmax - 1))[maskNjk], prob = par[skipindex], log = TRUE)
              Fv[j, maskNeik] <- pgeom(q = (0:(Kmax - 1))[maskNeik], prob = par[skipindex], lower.tail = FALSE, log.p = TRUE)
              skipindex <- skipindex + 1
            } else if (distr[j] == "nbinom") {
              fv[j, maskNjk] <- dnbinom(x = (0:(Kmax - 1))[maskNjk], size = par[skipindex], prob = par[skipindex + 1], log = TRUE)
              Fv[j, maskNeik] <- pnbinom(q = (0:(Kmax - 1))[maskNeik], size = par[skipindex], prob = par[skipindex + 1], lower.tail = FALSE, log.p = TRUE)
              skipindex <- skipindex + 2
            } else if (distr[j] == "pois") {
              fv[j, maskNjk] <- dpois(x = (0:(Kmax - 1))[maskNjk], lambda = par[skipindex], log = TRUE)
              Fv[j, maskNeik] <- ppois(q = (0:(Kmax - 1))[maskNeik], lambda = par[skipindex], lower.tail = FALSE, log.p = TRUE)
              skipindex <- skipindex + 1
            }
          }
          
          return(
            -(sum(res$Njk[1, ] * fv[1, ]) + sum(res$Njk[2, ] * fv[2, ])
              + sum(res$Neik[1, ] * Fv[2, ]) + sum(res$Neik[2, ] * Fv[1, ]))
          )
        }
        
      }
      
      # Constraints about the values of the parameters
      
      # Parameters values > 0
      u0 <- diag(x = 1, nrow = length(theta0), ncol = length(theta0))
      c0 <- rep(0, length(theta0))
      
      # Constraints on the values of the parameters
      
      u1 <- matrix(0, nrow = length(theta0), ncol = length(theta0))
      skipindex <- 1
      rowstoremove <- c()
      for (j in 1:S) {
        if (distr[j] == "dweibull") {
          u1[skipindex, skipindex] <- -1
          rowstoremove <- c(rowstoremove, skipindex + 1)
          skipindex <- skipindex + 2
        } else if (distr[j] == "geom") {
          u1[skipindex, skipindex] <- -1
          skipindex <- skipindex + 1
        } else if (distr[j] == "nbinom") {
          u1[skipindex + 1, skipindex + 1] <- -1
          rowstoremove <- c(rowstoremove, skipindex)
          skipindex <- skipindex + 2
        } else if (distr[j] == "pois") {
          rowstoremove <- c(rowstoremove, skipindex)
          skipindex <- skipindex + 1
        }
      }
      
      if (!is.null(rowstoremove)) {
        u1 <- u1[-rowstoremove, ]
      }
      
      if (!is.null(nrow(u1))) {
        c1 <- rep(-1, nrow(u1))
      } else {
        c1 <- c(-1)
      }
      
      if (length(u1) != 0) {
        u2 <- rbind(u0, u1)
        c2 <- c(c0, c1)
      } else {
        u2 <- u0
        c2 <- c0
      }
      
      CO2 <-
        constrOptim(
          theta = theta0,
          f = loglik,
          ui = u2,
          ci = c2,
          method = "Nelder-Mead"
        )
      
      skipindex <- 1
      for (j in 1:S) {
        if (distr[j] %in% c("dweibull", "nbinom")) {
          param[j, ] <- CO2$par[skipindex:(skipindex + 1)]
          skipindex <- skipindex + 2
        } else if (distr[j] %in% c("geom", "pois")) {
          param[j, 1] <- CO2$par[skipindex]
          skipindex <- skipindex + 1
        } else if (distr[j] == "unif") {
          param[j, ] <- .fit.param.fj.unif(res, j, Kmax)
        }
      }
      
    }
    
    ##########################################################
    # Case S > 2
    ##########################################################
  } else {
    
    if (cens.beg) {
      
      loglik <- function(par) {
        
        parpuv <- matrix(par[1:(S * (S - 2))], nrow = S, ncol = S - 2, byrow = T)
        parpuv <- cbind(parpuv, 1 - apply(parpuv, 1, sum))
        
        # Let's rebuild the transition matrix puv
        parp <- matrix(data = 0, nrow = S, ncol = S)
        parp[row(parp) != col(parp)] <- t(parpuv)
        parp <- t(parp)
        
        
        fv <- matrix(data = 0, nrow = S, ncol = Kmax)
        Fv <- matrix(data = 0, nrow = S, ncol = Kmax)
        Fv2 <- matrix(data = 0, nrow = S, ncol = Kmax)
        skipindex <- S * (S - 2) + 1
        for (j in 1:S) {
          
          maskNjk <- res$Njk[j, ] != 0
          maskNbjk <- res$Nbjk[j, ] != 0
          
          if (distr[j] == "dweibull") {
            fv[j, maskNjk] <- log(ddweibull(x = (1:Kmax)[maskNjk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
            Fv[j, maskNbjk] <- log(1 - pdweibull(x = (1:Kmax)[maskNbjk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE) + .Machine$double.xmin)
            Fv2[j, ] <- 1 - pdweibull(x = (1:Kmax), q = par[skipindex], beta = par[skipindex + 1], zero = FALSE)
            skipindex <- skipindex + 2
          } else if (distr[j] == "geom") {
            fv[j, maskNjk] <- dgeom(x = (0:(Kmax - 1))[maskNjk], prob = par[skipindex], log = TRUE)
            Fv[j, maskNbjk] <- pgeom(q = (0:(Kmax - 1))[maskNbjk], prob = par[skipindex], lower.tail = FALSE, log.p = TRUE)
            Fv2[j, ] <- pgeom(q = 0:(Kmax - 1), prob = par[skipindex], lower.tail = FALSE)
            skipindex <- skipindex + 1
          } else if (distr[j] == "nbinom") {
            fv[j, maskNjk] <- dnbinom(x = (0:(Kmax - 1))[maskNjk], size = par[skipindex], prob = par[skipindex + 1], log = TRUE)
            Fv[j, maskNbjk] <- pnbinom(q = (0:(Kmax - 1))[maskNbjk], size = par[skipindex], prob = par[skipindex + 1], lower.tail = FALSE, log.p = TRUE)
            Fv2[j, ] <- pnbinom(q = 0:(Kmax - 1), size = par[skipindex], prob = par[skipindex + 1], lower.tail = FALSE)
            skipindex <- skipindex + 2
          } else if (distr[j] == "pois") {
            fv[j, maskNjk] <- dpois(x = (0:(Kmax - 1))[maskNjk], lambda = par[skipindex], log = TRUE)
            Fv[j, maskNbjk] <- ppois(q = (0:(Kmax - 1))[maskNbjk], lambda = par[skipindex], lower.tail = FALSE, log.p = TRUE)
            Fv2[j, ] <- ppois(q = 0:(Kmax - 1), lambda = par[skipindex], lower.tail = FALSE)
            skipindex <- skipindex + 1
          }
        }
        
        Fv2 <- parp %*% Fv2
        Fv2[Fv2 != 0] <- log(Fv2[Fv2 != 0])
        
        return(-(
          sum(res$Nij[row(res$Nij) != col(res$Nij)] * log(parp[row(parp) != col(parp)])) +
            sum(res$Njk * fv) + sum(res$Nbjk * Fv) + sum(res$Neik * Fv2)
        ))
      }
      
    } else {
      
      loglik <- function(par) {
        
        parpuv <- matrix(par[1:(S * (S - 2))], nrow = S, ncol = S - 2, byrow = T)
        parpuv <- cbind(parpuv, 1 - apply(parpuv, 1, sum))
        
        # Let's rebuild the transition matrix puv
        parp <- matrix(data = 0, nrow = S, ncol = S)
        parp[row(parp) != col(parp)] <- t(parpuv)
        parp <- t(parp)
        
        
        fv <- matrix(data = 0, nrow = S, ncol = Kmax)
        Fv2 <- matrix(data = 0, nrow = S, ncol = Kmax)
        skipindex <- S * (S - 2) + 1
        for (j in 1:S) {
          
          maskNjk <- res$Njk[j, ] != 0
          
          if (distr[j] == "dweibull") {
            fv[j, maskNjk] <- log(ddweibull(x = (1:Kmax)[maskNjk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
            Fv2[j, ] <- 1 - pdweibull(x = 1:Kmax, q = par[skipindex], beta = par[skipindex + 1], zero = FALSE)
            skipindex <- skipindex + 2
          } else if (distr[j] == "geom") {
            fv[j, maskNjk] <- dgeom(x = (0:(Kmax - 1))[maskNjk], prob = par[skipindex], log = TRUE)
            Fv2[j, ] <- pgeom(q = 0:(Kmax - 1), prob = par[skipindex], lower.tail = FALSE)
            skipindex <- skipindex + 1
          } else if (distr[j] == "nbinom") {
            fv[j, maskNjk] <- dnbinom(x = (0:(Kmax - 1))[maskNjk], size = par[skipindex], prob = par[skipindex + 1], log = TRUE)
            Fv2[j, ] <- pnbinom(q = 0:(Kmax - 1), size = par[skipindex], prob = par[skipindex + 1], lower.tail = FALSE)
            skipindex <- skipindex + 2
          } else if (distr[j] == "pois") {
            fv[j, maskNjk] <- dpois(x = (0:(Kmax - 1))[maskNjk], lambda = par[skipindex], log = TRUE)
            Fv2[j, ] <- ppois(q = 0:(Kmax - 1), lambda = par[skipindex], lower.tail = FALSE)
            skipindex <- skipindex + 1
          }
        }
        
        Fv2 <- parp %*% Fv2
        Fv2[Fv2 != 0] <- log(Fv2[Fv2 != 0])
        
        return(-(
          sum(res$Nij[row(res$Nij) != col(res$Nij)] * log(parp[row(parp) != col(parp)])) +
            sum(res$Njk * fv) + sum(res$Neik * Fv2)
        ))
      }
      
    }
    
    # Constraints about the values of the parameters:
    
    # puv >= 0 and parameters values > 0
    u0 <- diag(x = 1, nrow = S * (S - 2) + length(theta0), ncol = S * (S - 2) + length(theta0))
    c0 <- rep(0, S * (S - 2) + length(theta0))
    
    # puv <= 1
    u1 <- matrix(0, nrow = S * (S - 2), ncol = S * (S - 2) + length(theta0))
    diag(u1) <- -1
    c1 <- rep(-1, S * (S - 2))
    
    # sum(puv) <= 1
    u2 <- matrix(0, nrow = S, ncol = S * (S - 2) + length(theta0))
    u2[1, 1:(S - 2)] <- -1
    for (l in 1:(S - 1)) {
      u2[l + 1, (l * (S - 2) + 1):(l * (S - 2) + (S - 2))] <- -1
    }
    c2 <- rep(-1, S)
    
    # Specific constraints depending on the distribution
    u3 <- matrix(0, nrow = length(theta0), ncol = S * (S - 2) + length(theta0))
    skipindex <- 1
    rowstoremove <- c()
    for (j in 1:S) {
      if (distr[j] == "dweibull") {
        u3[skipindex, skipindex + S * (S - 2)] <- -1
        rowstoremove <- c(rowstoremove, skipindex + 1)
        skipindex <- skipindex + 2
      } else if (distr[j] == "geom") {
        u3[skipindex, skipindex + S * (S - 2)] <- -1
        skipindex <- skipindex + 1
      } else if (distr[j] == "nbinom") {
        u3[skipindex + 1, skipindex + 1 + S * (S - 2)] <- -1
        rowstoremove <- c(rowstoremove, skipindex)
        skipindex <- skipindex + 2
      } else if (distr[j] == "pois") {
        rowstoremove <- c(rowstoremove, skipindex)
        skipindex <- skipindex + 1
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
        theta = c(rep(1 / (S - 1), S * (S - 2)), theta0),
        f = loglik,
        ui = u4,
        ci = c4,
        method = "Nelder-Mead"
      )
    
    # Rebuild parameters
    parpuv <- matrix(CO2$par[1:(S * (S - 2))], nrow = S, ncol = S - 2, byrow = T)
    parpuv <- cbind(parpuv, 1 - apply(parpuv, 1, sum))
    
    ptrans <- matrix(data = 0, nrow = S, ncol = S)
    ptrans[row(ptrans) != col(ptrans)] <- t(parpuv)
    ptrans <- t(ptrans)
    
    param <- matrix(data = NA, nrow = S, ncol = 2)
    skipindex <- S * (S - 2) + 1
    for (j in 1:S) {
      if (distr[j] %in% c("dweibull", "nbinom")) {
        param[j, ] <- CO2$par[skipindex:(skipindex + 1)]
        skipindex <- skipindex + 2
      } else if (distr[j] %in% c("geom", "pois")) {
        param[j, 1] <- CO2$par[skipindex]
        skipindex <- skipindex + 1
      } else if (distr[j] == "unif") {
        param[j, ] <- .fit.param.fj.unif(res, j, Kmax)
      }
    }
    
  }
  
  return(list(ptrans = ptrans, param = param))
}