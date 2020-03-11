.fit.param <- function(seq, E, type.sojourn, distr, cens.end, cens.beg) {
  
  S <- length(E)
  
  nbseq <- length(seq)
  Kmax <- 0
  J <- list()
  T <- list()
  L <- list()
  
  for (m in 1:nbseq) {
    
    processes <- .getprocesses(seq[[m]], E)
    
    J[[m]] <- processes$J
    T[[m]] <- processes$T
    L[[m]] <- T[[m]] - c(1, T[[m]][-length(T[[m]] - 1)]) # Sojourn time
    Kmax <- max(Kmax, max(L[[m]])) # Maximal sojourn time
  }
  
  res <- .count(J, L, S, Kmax)
  
  if (type.sojourn == "fij") {
    
    if (cens.end) {# We can't decompose the loglikelihood as a sum of optimization problems
      
      ptrans <- matrix(0, nrow = S, ncol = S)
      param <- array(data = NA, dim = c(S, S, 2))
      
      for (i in 1:S) {
        estparam <- .fit.param.fij.endcensoring(res, i, S, Kmax, distr, cens.beg)
        ptrans[i, ] <- estparam$ptrans
        param[i, , ] <- estparam$param
      }
      
    } else {
      
      # Estimation of the transition matrix - component #1
      ptrans <- res$Nij / res$Ni
      ptrans[which(is.na(ptrans))] <- 0
      
      param <- array(data = NA, dim = c(S, S, 2))
      
      for (i in 1:S) {
        for (j in 1:S) {
          if (i != j) {
            if (distr[i, j] == "dweibull") {
              param[i, j, ] <- .fit.param.fij.dweibull(res, i, j, Kmax, cens.beg)
            } else if (distr[i, j] == "geom") {
              param[i, j, ] <- .fit.param.fij.geom(res, i, j, Kmax, cens.beg)
            } else if (distr[i, j] == "nbinom") {
              param[i, j, ] <- .fit.param.fij.nbinom(res, i, j, Kmax, cens.beg)
            } else if (distr[i, j] == "pois") {
              param[i, j, ] <- .fit.param.fij.pois(res, i, j, Kmax, cens.beg)
            } else if (distr[i, j] == "unif") {
              param[i, j, ] <- .fit.param.fij.unif(res, i, j, Kmax)
            }
          }
        }
      }
      
    }
    
  } else if (type.sojourn == "fi") {
    
    # Estimation of the transition matrix - component #1
    ptrans <- res$Nij / res$Ni
    ptrans[which(is.na(ptrans))] <- 0
    
    param <- matrix(data = NA, nrow = S, ncol = 2)
    
    for (i in 1:S) {
      
      if (distr[i] == "dweibull") {
        param[i, ] <- .fit.param.fi.dweibull(res, i, Kmax, cens.beg, cens.end)
      } else if (distr[i] == "geom") {
        param[i, ] <- .fit.param.fi.geom(res, i, Kmax, cens.beg, cens.end)
      } else if (distr[i] == "nbinom") {
        param[i, ] <- .fit.param.fi.nbinom(res, i, Kmax, cens.beg, cens.end)
      } else if (distr[i] == "pois") {
        param[i, ] <- .fit.param.fi.pois(res, i, Kmax, cens.beg, cens.end)
      } else if (distr[i] == "unif") {
        param[i, ] <- .fit.param.fi.unif(res, i, Kmax)
      }
      
    }
    
  } else if (type.sojourn == "fj") {
    
    if (cens.end) {# We can't decompose the loglikelihood as a sum of optimization problems
      
      estparam <- .fit.param.fj.endcensoring(res, S, Kmax, distr, cens.beg)
      ptrans <- estparam$ptrans
      param <- estparam$param
      
    } else {
      
      # Estimation of the transition matrix - component #1
      ptrans <- res$Nij / res$Ni
      ptrans[which(is.na(ptrans))] <- 0
      
      param <- matrix(data = NA, nrow = S, ncol = 2)
      
      for (j in 1:S) {
        
        if (distr[j] == "dweibull") {
          param[j, ] <- .fit.param.fj.dweibull(res, j, Kmax, cens.beg)
        } else if (distr[j] == "geom") {
          param[j, ] <- .fit.param.fj.geom(res, j, Kmax, cens.beg)
        } else if (distr[j] == "nbinom") {
          param[j, ] <- .fit.param.fj.nbinom(res, j, Kmax, cens.beg)
        } else if (distr[j] == "pois") {
          param[j, ] <- .fit.param.fj.pois(res, j, Kmax, cens.beg)
        } else if (distr[j] == "unif") {
          param[j, ] <- .fit.param.fj.unif(res, j, Kmax)
        }
        
      }
      
    }
    
  } else if (type.sojourn == "f") {
    
    # Estimation of the transition matrix
    ptrans <- res$Nij / res$Ni
    ptrans[which(is.na(ptrans))] <- 0
    
    if (distr == "dweibull") {
      param <- .fit.param.f.dweibull(res, Kmax, cens.beg, cens.end)
    } else if (distr == "geom") {
      param <- .fit.param.f.geom(res, Kmax, cens.beg, cens.end)
    } else if (distr == "nbinom") {
      param <- .fit.param.f.nbinom(res, Kmax, cens.beg, cens.end)
    } else if (distr == "pois") {
      param <- .fit.param.f.pois(res, Kmax, cens.beg, cens.end)
    } else if (distr == "unif") {
      param <- .fit.param.f.unif(res, Kmax)
    }
    
  }
  
  # Inital distribution
  if (nbseq >= S * 10) {
    init <- res$Nstarti / sum(res$Nstarti)
  } else {# Computation of the limit distribution
    q <- .get.kernel(type.sojourn = type.sojourn, distr = distr, param = param, pij = ptrans, Kmax = Kmax, S = S)
    init <- .limit.distribution(q = q, pij = ptrans)
  }
  
  estimate <-
    list(
      E = E,
      init = init,
      type.sojourn = type.sojourn,
      ptrans = ptrans,
      distr = distr,
      param = param,
      cens.beg = FALSE,
      cens.end = FALSE
    )
  
  return(estimate)
  
}
