.fit.param <- function(sequences, states, type.sojourn, distr, cens.end, cens.beg) {
  
  states <- sequences$states
  s <- sequences$s
  L <- sequences$L
  kmax <- sequences$kmax
  counting <- sequences$counting
  
  
  if (type.sojourn == "fij") {
    
    if (cens.end) {# We can't decompose the loglikelihood as a sum of optimization problems
      
      ptrans <- matrix(0, nrow = s, ncol = s)
      param <- array(data = NA, dim = c(s, s, 2))
      
      for (i in 1:s) {
        estparam <- .fit.param.fij.endcensoring(counting, i, s, kmax, distr, cens.beg)
        ptrans[i, ] <- estparam$ptrans
        param[i, , ] <- estparam$param
      }
      
    } else {
      
      # Estimation of the transition matrix - component #1
      ptrans <- counting$Nij / counting$Ni
      ptrans[which(is.na(ptrans))] <- 0
      
      # Estimation of the parameters of each distribution
      param <- array(data = NA, dim = c(s, s, 2))
      
      for (i in 1:s) {
        for (j in 1:s) {
          if (i != j & !is.na(distr[i, j])) {
            if (distr[i, j] == "dweibull") {
              param[i, j, ] <- .fit.param.fij.dweibull(counting, i, j, kmax, cens.beg)
            } else if (distr[i, j] == "geom") {
              param[i, j, ] <- .fit.param.fij.geom(counting, i, j, kmax, cens.beg)
            } else if (distr[i, j] == "nbinom") {
              param[i, j, ] <- .fit.param.fij.nbinom(counting, i, j, kmax, cens.beg)
            } else if (distr[i, j] == "pois") {
              param[i, j, ] <- .fit.param.fij.pois(counting, i, j, kmax, cens.beg)
            } else if (distr[i, j] == "unif") {
              param[i, j, ] <- .fit.param.fij.unif(counting, i, j, kmax)
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
    param <- matrix(data = NA, nrow = s, ncol = 2)
    
    for (i in 1:s) {
      if (!is.na(distr[i])) {
        if (distr[i] == "dweibull") {
          param[i, ] <- .fit.param.fi.dweibull(counting, i, kmax, cens.beg, cens.end)
        } else if (distr[i] == "geom") {
          param[i, ] <- .fit.param.fi.geom(counting, i, kmax, cens.beg, cens.end)
        } else if (distr[i] == "nbinom") {
          param[i, ] <- .fit.param.fi.nbinom(counting, i, kmax, cens.beg, cens.end)
        } else if (distr[i] == "pois") {
          param[i, ] <- .fit.param.fi.pois(counting, i, kmax, cens.beg, cens.end)
        } else if (distr[i] == "unif") {
          param[i, ] <- .fit.param.fi.unif(counting, i, kmax)
        }
      }
    }
    
  } else if (type.sojourn == "fj") {
    
    if (cens.end) {# We can't decompose the loglikelihood as a sum of optimization problems
      
      estparam <- .fit.param.fj.endcensoring(counting, s, kmax, distr, cens.beg)
      ptrans <- estparam$ptrans
      param <- estparam$param
      
    } else {
      
      # Estimation of the transition matrix - component #1
      ptrans <- counting$Nij / counting$Ni
      ptrans[which(is.na(ptrans))] <- 0
      
      # Estimation of the parameters of each distribution
      param <- matrix(data = NA, nrow = s, ncol = 2)
      
      for (j in 1:s) {
        if (!is.na(distr[j])) {
          if (distr[j] == "dweibull") {
            param[j, ] <- .fit.param.fj.dweibull(counting, j, kmax, cens.beg)
          } else if (distr[j] == "geom") {
            param[j, ] <- .fit.param.fj.geom(counting, j, kmax, cens.beg)
          } else if (distr[j] == "nbinom") {
            param[j, ] <- .fit.param.fj.nbinom(counting, j, kmax, cens.beg)
          } else if (distr[j] == "pois") {
            param[j, ] <- .fit.param.fj.pois(counting, j, kmax, cens.beg)
          } else if (distr[j] == "unif") {
            param[j, ] <- .fit.param.fj.unif(counting, j, kmax)
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
      param <- .fit.param.f.dweibull(counting, kmax, cens.beg, cens.end)
    } else if (distr == "geom") {
      param <- .fit.param.f.geom(counting, kmax, cens.beg, cens.end)
    } else if (distr == "nbinom") {
      param <- .fit.param.f.nbinom(counting, kmax, cens.beg, cens.end)
    } else if (distr == "pois") {
      param <- .fit.param.f.pois(counting, kmax, cens.beg, cens.end)
    } else if (distr == "unif") {
      param <- .fit.param.f.unif(counting, kmax)
    }
    
  }
  
  estimate <-
    list(
      states = states,
      s = s,
      init = rep.int(x = NA, times = s),
      type.sojourn = type.sojourn,
      ptrans = ptrans,
      distr = distr,
      param = param,
      cens.beg = cens.beg,
      cens.end = cens.end
    )
  
  class(estimate) <- c("smm", "smmparametric")
  
  # Inital distribution
  if (L >= s * 10) {
    estimate$init <- counting$Nstarti / sum(counting$Nstarti)
  } else {# Computation of the limit distribution
    q <- .get.q(estimate, kmax)
    estimate$init <- .limitDistribution(q = q, ptrans = estimate$ptrans)
  }
  
  return(estimate)
  
}
