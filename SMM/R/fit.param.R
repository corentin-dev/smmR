.fit.param <- function(seq, E, type.sojourn, distr, cens.end, cens.beg) {
  
  S <- seq$S
  E <- seq$E
  L <- seq$L
  S <- seq$S
  Kmax <- seq$Kmax
  counting <- seq$counting
  
  
  if (type.sojourn == "fij") {
    
    if (cens.end) {# We can't decompose the loglikelihood as a sum of optimization problems
      
      ptrans <- matrix(0, nrow = S, ncol = S)
      param <- array(data = NA, dim = c(S, S, 2))
      
      for (i in 1:S) {
        estparam <- .fit.param.fij.endcensoring(counting, i, S, Kmax, distr, cens.beg)
        ptrans[i, ] <- estparam$ptrans
        param[i, , ] <- estparam$param
      }
      
    } else {
      
      # Estimation of the transition matrix - component #1
      ptrans <- counting$Nij / counting$Ni
      ptrans[which(is.na(ptrans))] <- 0
      
      # Estimation of the parameters of each distribution
      param <- array(data = NA, dim = c(S, S, 2))
      
      for (i in 1:S) {
        for (j in 1:S) {
          if (i != j & !is.na(distr[i, j])) {
            if (distr[i, j] == "dweibull") {
              param[i, j, ] <- .fit.param.fij.dweibull(counting, i, j, Kmax, cens.beg)
            } else if (distr[i, j] == "geom") {
              param[i, j, ] <- .fit.param.fij.geom(counting, i, j, Kmax, cens.beg)
            } else if (distr[i, j] == "nbinom") {
              param[i, j, ] <- .fit.param.fij.nbinom(counting, i, j, Kmax, cens.beg)
            } else if (distr[i, j] == "pois") {
              param[i, j, ] <- .fit.param.fij.pois(counting, i, j, Kmax, cens.beg)
            } else if (distr[i, j] == "unif") {
              param[i, j, ] <- .fit.param.fij.unif(counting, i, j, Kmax)
            }
          }
        }
      }
      
    }
    
  } else if (type.sojourn == "fi") {
    
    # Estimation of the transition matrix - component #1
    ptrans <- counting$Nij / counting$Ni
    ptrans[which(is.na(ptrans))] <- 0
    
    # Estimation of the parameters of each distribution
    param <- matrix(data = NA, nrow = S, ncol = 2)
    
    for (i in 1:S) {
      if (!is.na(distr[i])) {
        if (distr[i] == "dweibull") {
          param[i, ] <- .fit.param.fi.dweibull(counting, i, Kmax, cens.beg, cens.end)
        } else if (distr[i] == "geom") {
          param[i, ] <- .fit.param.fi.geom(counting, i, Kmax, cens.beg, cens.end)
        } else if (distr[i] == "nbinom") {
          param[i, ] <- .fit.param.fi.nbinom(counting, i, Kmax, cens.beg, cens.end)
        } else if (distr[i] == "pois") {
          param[i, ] <- .fit.param.fi.pois(counting, i, Kmax, cens.beg, cens.end)
        } else if (distr[i] == "unif") {
          param[i, ] <- .fit.param.fi.unif(counting, i, Kmax)
        }
      }
    }
    
  } else if (type.sojourn == "fj") {
    
    if (cens.end) {# We can't decompose the loglikelihood as a sum of optimization problems
      
      estparam <- .fit.param.fj.endcensoring(counting, S, Kmax, distr, cens.beg)
      ptrans <- estparam$ptrans
      param <- estparam$param
      
    } else {
      
      # Estimation of the transition matrix - component #1
      ptrans <- counting$Nij / counting$Ni
      ptrans[which(is.na(ptrans))] <- 0
      
      # Estimation of the parameters of each distribution
      param <- matrix(data = NA, nrow = S, ncol = 2)
      
      for (j in 1:S) {
        if (!is.na(distr[j])) {
          if (distr[j] == "dweibull") {
            param[j, ] <- .fit.param.fj.dweibull(counting, j, Kmax, cens.beg)
          } else if (distr[j] == "geom") {
            param[j, ] <- .fit.param.fj.geom(counting, j, Kmax, cens.beg)
          } else if (distr[j] == "nbinom") {
            param[j, ] <- .fit.param.fj.nbinom(counting, j, Kmax, cens.beg)
          } else if (distr[j] == "pois") {
            param[j, ] <- .fit.param.fj.pois(counting, j, Kmax, cens.beg)
          } else if (distr[j] == "unif") {
            param[j, ] <- .fit.param.fj.unif(counting, j, Kmax)
          }
        }
      }
      
    }
    
  } else if (type.sojourn == "f") {
    
    # Estimation of the transition matrix
    ptrans <- counting$Nij / counting$Ni
    ptrans[which(is.na(ptrans))] <- 0
    
    # Estimation of the parameters of each distribution
    param <- rep.int(x = NA, times = 2)
    
    if (distr == "dweibull") {
      param <- .fit.param.f.dweibull(counting, Kmax, cens.beg, cens.end)
    } else if (distr == "geom") {
      param <- .fit.param.f.geom(counting, Kmax, cens.beg, cens.end)
    } else if (distr == "nbinom") {
      param <- .fit.param.f.nbinom(counting, Kmax, cens.beg, cens.end)
    } else if (distr == "pois") {
      param <- .fit.param.f.pois(counting, Kmax, cens.beg, cens.end)
    } else if (distr == "unif") {
      param <- .fit.param.f.unif(counting, Kmax)
    }
    
  }
  
  estimate <-
    list(
      E = E,
      S = S,
      init = rep.int(x = NA, times = S),
      type.sojourn = type.sojourn,
      ptrans = ptrans,
      distr = distr,
      param = param,
      cens.beg = cens.beg,
      cens.end = cens.end
    )
  
  class(estimate) <- c("smm", "smmparametric")
  
  # Inital distribution
  if (L >= S * 10) {
    estimate$init <- counting$Nstarti / sum(counting$Nstarti)
  } else {# Computation of the limit distribution
    q <- .get.q(estimate, Kmax)
    estimate$init <- .limitDistribution(q = q, ptrans = estimate$ptrans)
  }
  
  return(estimate)
  
}
