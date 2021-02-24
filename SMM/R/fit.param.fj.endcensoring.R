.fit.param.fj.endcensoring <- function(counting, s, kmax, distr, cens.beg) {
  
  # We use the non-censoring case as initial values for the optimization problem
  theta0 <- matrix(data = NA, nrow = 2, ncol = s)
  for (j in 1:s) {
    
    if (distr[j] == "dweibull") {
      theta0[, j] <- .fit.param.fj.dweibull(counting, j, kmax, cens.beg = FALSE)
    } else if (distr[j] == "geom") {
      theta0[, j] <- .fit.param.fj.geom(counting, j, kmax, cens.beg = FALSE)
    } else if (distr[j] == "nbinom") {
      theta0[, j] <- .fit.param.fj.nbinom(counting, j, kmax, cens.beg = FALSE)
    } else if (distr[j] == "pois") {
      theta0[, j] <- .fit.param.fj.pois(counting, j, kmax, cens.beg = FALSE)
    }
    
  }
  
  theta0 <- as.vector(theta0[!(is.na(theta0))])
  
  ##########################################################
  # Specific case s = 2
  ##########################################################
  if (s == 2) {
    
    ptrans <- matrix(c(0, 1, 1, 0), nrow = s, byrow = TRUE)
    param <- matrix(data = NA, nrow = s, ncol = 2)
    
    ##########################################################
    # Sojourn time distributions are 2 uniforms
    ##########################################################
    if ((distr[1] == "unif") & (distr[2] == "unif")) {
      
      for (j in 1:s) {
        param[j, ] <- .fit.param.fj.unif(counting, j, kmax)
      }
      
      ##########################################################
      # There is one uniform distribution among the 2 distributions
      ##########################################################
    } else if ("unif" %in% distr) {
      
      for (j in 1:s) {
        if (distr[j] == "unif") {
          param[j, ] <- .fit.param.fj.unif(counting, j, kmax)
        } else {
          
          if (cens.beg) {
            
            loglik <- function(par) {
              
              fv <- matrix(data = 0, nrow = s, ncol = kmax)
              Fv <- matrix(data = 0, nrow = s, ncol = kmax)
              skipindex <- 1
              for (j in 1:s) {
                
                maskNjk <- counting$Njk[j, ] != 0
                maskNeNb <- (counting$Nbjk[j, ] + counting$Neik[abs(j - 3), ]) != 0
                
                if (distr[j] == "dweibull") {
                  fv[j, maskNjk] <- log(ddweibull(x = (1:kmax)[maskNjk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
                  Fv[j, maskNeNb] <- log(1 - pdweibull(x = (1:kmax)[maskNeNb], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE) + .Machine$double.xmin)
                  skipindex <- skipindex + 2
                } else if (distr[j] == "geom") {
                  fv[j, maskNjk] <- dgeom(x = (0:(kmax - 1))[maskNjk], prob = par[skipindex], log = TRUE)
                  Fv[j, maskNeNb] <- pgeom(q = (0:(kmax - 1))[maskNeNb], prob = par[skipindex], lower.tail = FALSE, log.p = TRUE)
                  skipindex <- skipindex + 1
                } else if (distr[j] == "nbinom") {
                  fv[j, maskNjk] <- dnbinom(x = (0:(kmax - 1))[maskNjk], size = par[skipindex], prob = par[skipindex + 1], log = TRUE)
                  Fv[j, maskNeNb] <- pnbinom(q = (0:(kmax - 1))[maskNeNb], size = par[skipindex], prob = par[skipindex + 1], lower.tail = FALSE, log.p = TRUE)
                  skipindex <- skipindex + 2
                } else if (distr[j] == "pois") {
                  fv[j, maskNjk] <- dpois(x = (0:(kmax - 1))[maskNjk], lambda = par[skipindex], log = TRUE)
                  Fv[j, maskNeNb] <- ppois(q = (0:(kmax - 1))[maskNeNb], lambda = par[skipindex], lower.tail = FALSE, log.p = TRUE)
                  skipindex <- skipindex + 1
                }
              }
              
              return(
                -(sum(counting$Njk[1, ] * fv[1, ]) + sum(counting$Njk[2, ] * fv[2, ])
                  + sum((counting$Nbjk[1, ] + counting$Neik[2, ]) * Fv[1, ]) 
                  + sum((counting$Nbjk[2, ] + counting$Neik[1, ]) * Fv[2, ]))
              )
            }
            
          } else {
            
            loglik <- function(par) {
              
              fv <- matrix(data = 0, nrow = s, ncol = kmax)
              Fv <- matrix(data = 0, nrow = s, ncol = kmax)
              skipindex <- 1
              for (j in 1:s) {
                
                maskNjk <- counting$Njk[j, ] != 0
                maskNeik <- counting$Neik[abs(j - 3), ] != 0
                
                if (distr[j] == "dweibull") {
                  fv[j, maskNjk] <- log(ddweibull(x = (1:kmax)[maskNjk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
                  Fv[j, maskNeik] <- log(1 - pdweibull(x = (1:kmax)[maskNeik], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE) + .Machine$double.xmin)
                  skipindex <- skipindex + 2
                } else if (distr[j] == "geom") {
                  fv[j, maskNjk] <- dgeom(x = (0:(kmax - 1))[maskNjk], prob = par[skipindex], log = TRUE)
                  Fv[j, maskNeik] <- pgeom(q = (0:(kmax - 1))[maskNeik], prob = par[skipindex], lower.tail = FALSE, log.p = TRUE)
                  skipindex <- skipindex + 1
                } else if (distr[j] == "nbinom") {
                  fv[j, maskNjk] <- dnbinom(x = (0:(kmax - 1))[maskNjk], size = par[skipindex], prob = par[skipindex + 1], log = TRUE)
                  Fv[j, maskNeik] <- pnbinom(q = (0:(kmax - 1))[maskNeik], size = par[skipindex], prob = par[skipindex + 1], lower.tail = FALSE, log.p = TRUE)
                  skipindex <- skipindex + 2
                } else if (distr[j] == "pois") {
                  fv[j, maskNjk] <- dpois(x = (0:(kmax - 1))[maskNjk], lambda = par[skipindex], log = TRUE)
                  Fv[j, maskNeik] <- ppois(q = (0:(kmax - 1))[maskNeik], lambda = par[skipindex], lower.tail = FALSE, log.p = TRUE)
                  skipindex <- skipindex + 1
                }
              }
              
              return(
                -(sum(counting$Njk[1, ] * fv[1, ]) + sum(counting$Njk[2, ] * fv[2, ])
                  + sum(counting$Neik[1, ] * Fv[2, ]) + sum(counting$Neik[2, ] * Fv[1, ]))
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
            
            mle <- constrOptim(
              theta = theta0,
              f = loglik,
              ui = rbind(u0, u1),
              ci = c(c0, c1),
              method = "Nelder-Mead"
            )
            
            param[j, ] <- mle$par
            
          } else if (distr[j] == "geom") {
            
            mle <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = 1)
            param[j, ] <- mle$par
            
          } else if (distr[j] == "nbinom") {
            
            # Constraints about the values of the parameters
            
            # alpha, p > 0
            u0 <- diag(x = 1, nrow = 2)
            c0 <- c(0, 0)
            
            # p < 1
            u1 <- matrix(data = c(0, -1), nrow = 1, ncol = 2)
            c1 <- c(-1)
            
            mle <- constrOptim(
              theta = theta0,
              f = loglik,
              ui = rbind(u0, u1),
              ci = c(c0, c1),
              method = "Nelder-Mead"
            )
            param[j, ] <- mle$par
            
          } else if (distr[j] == "pois") {
            
            mle <- optim(par = theta0, loglik, method = "Brent", lower = 0, upper = kmax - 1)
            param[j, ] <- mle$par
            
          }
        }
      }
      
      ##########################################################
      # Two sojourn time distributions different from the uniform
      ##########################################################
    } else {
      
      if (cens.beg) {
        
        loglik <- function(par) {
          
          fv <- matrix(data = 0, nrow = s, ncol = kmax)
          Fv <- matrix(data = 0, nrow = s, ncol = kmax)
          skipindex <- 1
          for (j in 1:s) {
            
            maskNjk <- counting$Njk[j, ] != 0
            maskNeNb <- (counting$Nbjk[j, ] + counting$Neik[abs(j - 3), ]) != 0
            
            if (distr[j] == "dweibull") {
              fv[j, maskNjk] <- log(ddweibull(x = (1:kmax)[maskNjk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
              Fv[j, maskNeNb] <- log(1 - pdweibull(x = (1:kmax)[maskNeNb], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE) + .Machine$double.xmin)
              skipindex <- skipindex + 2
            } else if (distr[j] == "geom") {
              fv[j, maskNjk] <- dgeom(x = (0:(kmax - 1))[maskNjk], prob = par[skipindex], log = TRUE)
              Fv[j, maskNeNb] <- pgeom(q = (0:(kmax - 1))[maskNeNb], prob = par[skipindex], lower.tail = FALSE, log.p = TRUE)
              skipindex <- skipindex + 1
            } else if (distr[j] == "nbinom") {
              fv[j, maskNjk] <- dnbinom(x = (0:(kmax - 1))[maskNjk], size = par[skipindex], prob = par[skipindex + 1], log = TRUE)
              Fv[j, maskNeNb] <- pnbinom(q = (0:(kmax - 1))[maskNeNb], size = par[skipindex], prob = par[skipindex + 1], lower.tail = FALSE, log.p = TRUE)
              skipindex <- skipindex + 2
            } else if (distr[j] == "pois") {
              fv[j, maskNjk] <- dpois(x = (0:(kmax - 1))[maskNjk], lambda = par[skipindex], log = TRUE)
              Fv[j, maskNeNb] <- ppois(q = (0:(kmax - 1))[maskNeNb], lambda = par[skipindex], lower.tail = FALSE, log.p = TRUE)
              skipindex <- skipindex + 1
            }
          }
          
          return(
            -(sum(counting$Njk[1, ] * fv[1, ]) + sum(counting$Njk[2, ] * fv[2, ])
              + sum((counting$Nbjk[1, ] + counting$Neik[2, ]) * Fv[1, ]) 
              + sum((counting$Nbjk[2, ] + counting$Neik[1, ]) * Fv[2, ])
            )
          )
        }
        
      } else {
        
        loglik <- function(par) {
          
          fv <- matrix(data = 0, nrow = s, ncol = kmax)
          Fv <- matrix(data = 0, nrow = s, ncol = kmax)
          skipindex <- 1
          for (j in 1:s) {
            
            maskNjk <- counting$Njk[j, ] != 0
            maskNeik <- counting$Neik[abs(j - 3), ] != 0
            
            if (distr[j] == "dweibull") {
              fv[j, maskNjk] <- log(ddweibull(x = (1:kmax)[maskNjk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
              Fv[j, maskNeik] <- log(1 - pdweibull(x = (1:kmax)[maskNeik], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE) + .Machine$double.xmin)
              skipindex <- skipindex + 2
            } else if (distr[j] == "geom") {
              fv[j, maskNjk] <- dgeom(x = (0:(kmax - 1))[maskNjk], prob = par[skipindex], log = TRUE)
              Fv[j, maskNeik] <- pgeom(q = (0:(kmax - 1))[maskNeik], prob = par[skipindex], lower.tail = FALSE, log.p = TRUE)
              skipindex <- skipindex + 1
            } else if (distr[j] == "nbinom") {
              fv[j, maskNjk] <- dnbinom(x = (0:(kmax - 1))[maskNjk], size = par[skipindex], prob = par[skipindex + 1], log = TRUE)
              Fv[j, maskNeik] <- pnbinom(q = (0:(kmax - 1))[maskNeik], size = par[skipindex], prob = par[skipindex + 1], lower.tail = FALSE, log.p = TRUE)
              skipindex <- skipindex + 2
            } else if (distr[j] == "pois") {
              fv[j, maskNjk] <- dpois(x = (0:(kmax - 1))[maskNjk], lambda = par[skipindex], log = TRUE)
              Fv[j, maskNeik] <- ppois(q = (0:(kmax - 1))[maskNeik], lambda = par[skipindex], lower.tail = FALSE, log.p = TRUE)
              skipindex <- skipindex + 1
            }
          }
          
          return(
            -(sum(counting$Njk[1, ] * fv[1, ]) + sum(counting$Njk[2, ] * fv[2, ])
              + sum(counting$Neik[1, ] * Fv[2, ]) + sum(counting$Neik[2, ] * Fv[1, ]))
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
      for (j in 1:s) {
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
      
      mle <-
        constrOptim(
          theta = theta0,
          f = loglik,
          ui = u2,
          ci = c2,
          method = "Nelder-Mead"
        )
      
      skipindex <- 1
      for (j in 1:s) {
        if (distr[j] %in% c("dweibull", "nbinom")) {
          param[j, ] <- mle$par[skipindex:(skipindex + 1)]
          skipindex <- skipindex + 2
        } else if (distr[j] %in% c("geom", "pois")) {
          param[j, 1] <- mle$par[skipindex]
          skipindex <- skipindex + 1
        } else if (distr[j] == "unif") {
          param[j, ] <- .fit.param.fj.unif(counting, j, kmax)
        }
      }
      
    }
    
    ##########################################################
    # Case s > 2
    ##########################################################
  } else {
    
    if (cens.beg) {
      
      loglik <- function(par) {
        
        parpuv <- matrix(par[1:(s * (s - 2))], nrow = s, ncol = s - 2, byrow = T)
        parpuv <- cbind(parpuv, 1 - apply(parpuv, 1, sum))
        
        # Let's rebuild the transition matrix puv
        parp <- matrix(data = 0, nrow = s, ncol = s)
        parp[row(parp) != col(parp)] <- t(parpuv)
        parp <- t(parp)
        
        
        fv <- matrix(data = 0, nrow = s, ncol = kmax)
        Fv <- matrix(data = 0, nrow = s, ncol = kmax)
        Fv2 <- matrix(data = 0, nrow = s, ncol = kmax)
        skipindex <- s * (s - 2) + 1
        for (j in 1:s) {
          
          maskNjk <- counting$Njk[j, ] != 0
          maskNbjk <- counting$Nbjk[j, ] != 0
          
          if (distr[j] == "dweibull") {
            fv[j, maskNjk] <- log(ddweibull(x = (1:kmax)[maskNjk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
            Fv[j, maskNbjk] <- log(1 - pdweibull(x = (1:kmax)[maskNbjk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE) + .Machine$double.xmin)
            Fv2[j, ] <- 1 - pdweibull(x = (1:kmax), q = par[skipindex], beta = par[skipindex + 1], zero = FALSE)
            skipindex <- skipindex + 2
          } else if (distr[j] == "geom") {
            fv[j, maskNjk] <- dgeom(x = (0:(kmax - 1))[maskNjk], prob = par[skipindex], log = TRUE)
            Fv[j, maskNbjk] <- pgeom(q = (0:(kmax - 1))[maskNbjk], prob = par[skipindex], lower.tail = FALSE, log.p = TRUE)
            Fv2[j, ] <- pgeom(q = 0:(kmax - 1), prob = par[skipindex], lower.tail = FALSE)
            skipindex <- skipindex + 1
          } else if (distr[j] == "nbinom") {
            fv[j, maskNjk] <- dnbinom(x = (0:(kmax - 1))[maskNjk], size = par[skipindex], prob = par[skipindex + 1], log = TRUE)
            Fv[j, maskNbjk] <- pnbinom(q = (0:(kmax - 1))[maskNbjk], size = par[skipindex], prob = par[skipindex + 1], lower.tail = FALSE, log.p = TRUE)
            Fv2[j, ] <- pnbinom(q = 0:(kmax - 1), size = par[skipindex], prob = par[skipindex + 1], lower.tail = FALSE)
            skipindex <- skipindex + 2
          } else if (distr[j] == "pois") {
            fv[j, maskNjk] <- dpois(x = (0:(kmax - 1))[maskNjk], lambda = par[skipindex], log = TRUE)
            Fv[j, maskNbjk] <- ppois(q = (0:(kmax - 1))[maskNbjk], lambda = par[skipindex], lower.tail = FALSE, log.p = TRUE)
            Fv2[j, ] <- ppois(q = 0:(kmax - 1), lambda = par[skipindex], lower.tail = FALSE)
            skipindex <- skipindex + 1
          }
        }
        
        Fv2 <- parp %*% Fv2
        Fv2[Fv2 != 0] <- log(Fv2[Fv2 != 0])
        
        return(-(
          sum(counting$Nij[row(counting$Nij) != col(counting$Nij)] * log(parp[row(parp) != col(parp)])) +
            sum(counting$Njk * fv) + sum(counting$Nbjk * Fv) + sum(counting$Neik * Fv2)
        ))
      }
      
    } else {
      
      loglik <- function(par) {
        
        parpuv <- matrix(par[1:(s * (s - 2))], nrow = s, ncol = s - 2, byrow = T)
        parpuv <- cbind(parpuv, 1 - apply(parpuv, 1, sum))
        
        # Let's rebuild the transition matrix puv
        parp <- matrix(data = 0, nrow = s, ncol = s)
        parp[row(parp) != col(parp)] <- t(parpuv)
        parp <- t(parp)
        
        
        fv <- matrix(data = 0, nrow = s, ncol = kmax)
        Fv2 <- matrix(data = 0, nrow = s, ncol = kmax)
        skipindex <- s * (s - 2) + 1
        for (j in 1:s) {
          
          maskNjk <- counting$Njk[j, ] != 0
          
          if (distr[j] == "dweibull") {
            fv[j, maskNjk] <- log(ddweibull(x = (1:kmax)[maskNjk], q = par[skipindex], beta = par[skipindex + 1], zero = FALSE))
            Fv2[j, ] <- 1 - pdweibull(x = 1:kmax, q = par[skipindex], beta = par[skipindex + 1], zero = FALSE)
            skipindex <- skipindex + 2
          } else if (distr[j] == "geom") {
            fv[j, maskNjk] <- dgeom(x = (0:(kmax - 1))[maskNjk], prob = par[skipindex], log = TRUE)
            Fv2[j, ] <- pgeom(q = 0:(kmax - 1), prob = par[skipindex], lower.tail = FALSE)
            skipindex <- skipindex + 1
          } else if (distr[j] == "nbinom") {
            fv[j, maskNjk] <- dnbinom(x = (0:(kmax - 1))[maskNjk], size = par[skipindex], prob = par[skipindex + 1], log = TRUE)
            Fv2[j, ] <- pnbinom(q = 0:(kmax - 1), size = par[skipindex], prob = par[skipindex + 1], lower.tail = FALSE)
            skipindex <- skipindex + 2
          } else if (distr[j] == "pois") {
            fv[j, maskNjk] <- dpois(x = (0:(kmax - 1))[maskNjk], lambda = par[skipindex], log = TRUE)
            Fv2[j, ] <- ppois(q = 0:(kmax - 1), lambda = par[skipindex], lower.tail = FALSE)
            skipindex <- skipindex + 1
          }
        }
        
        Fv2 <- parp %*% Fv2
        Fv2[Fv2 != 0] <- log(Fv2[Fv2 != 0])
        
        return(-(
          sum(counting$Nij[row(counting$Nij) != col(counting$Nij)] * log(parp[row(parp) != col(parp)])) +
            sum(counting$Njk * fv) + sum(counting$Neik * Fv2)
        ))
      }
      
    }
    
    # Constraints about the values of the parameters:
    
    # puv >= 0 and parameters values > 0
    u0 <- diag(x = 1, nrow = s * (s - 2) + length(theta0), ncol = s * (s - 2) + length(theta0))
    c0 <- rep(0, s * (s - 2) + length(theta0))
    
    # puv <= 1
    u1 <- matrix(0, nrow = s * (s - 2), ncol = s * (s - 2) + length(theta0))
    diag(u1) <- -1
    c1 <- rep(-1, s * (s - 2))
    
    # sum(puv) <= 1
    u2 <- matrix(0, nrow = s, ncol = s * (s - 2) + length(theta0))
    u2[1, 1:(s - 2)] <- -1
    for (l in 1:(s - 1)) {
      u2[l + 1, (l * (s - 2) + 1):(l * (s - 2) + (s - 2))] <- -1
    }
    c2 <- rep(-1, s)
    
    # Specific constraints depending on the distribution
    u3 <- matrix(0, nrow = length(theta0), ncol = s * (s - 2) + length(theta0))
    skipindex <- 1
    rowstoremove <- c()
    for (j in 1:s) {
      if (distr[j] == "dweibull") {
        u3[skipindex, skipindex + s * (s - 2)] <- -1
        rowstoremove <- c(rowstoremove, skipindex + 1)
        skipindex <- skipindex + 2
      } else if (distr[j] == "geom") {
        u3[skipindex, skipindex + s * (s - 2)] <- -1
        skipindex <- skipindex + 1
      } else if (distr[j] == "nbinom") {
        u3[skipindex + 1, skipindex + 1 + s * (s - 2)] <- -1
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
    
    mle <-
      constrOptim(
        theta = c(rep(1 / (s - 1), s * (s - 2)), theta0),
        f = loglik,
        ui = u4,
        ci = c4,
        method = "Nelder-Mead"
      )
    
    # Rebuild parameters
    parpuv <- matrix(mle$par[1:(s * (s - 2))], nrow = s, ncol = s - 2, byrow = T)
    parpuv <- cbind(parpuv, 1 - apply(parpuv, 1, sum))
    
    ptrans <- matrix(data = 0, nrow = s, ncol = s)
    ptrans[row(ptrans) != col(ptrans)] <- t(parpuv)
    ptrans <- t(ptrans)
    
    param <- matrix(data = NA, nrow = s, ncol = 2)
    skipindex <- s * (s - 2) + 1
    for (j in 1:s) {
      if (distr[j] %in% c("dweibull", "nbinom")) {
        param[j, ] <- mle$par[skipindex:(skipindex + 1)]
        skipindex <- skipindex + 2
      } else if (distr[j] %in% c("geom", "pois")) {
        param[j, 1] <- mle$par[skipindex]
        skipindex <- skipindex + 1
      } else if (distr[j] == "unif") {
        param[j, ] <- .fit.param.fj.unif(counting, j, kmax)
      }
    }
    
  }
  
  return(list(ptrans = ptrans, param = param))
}
