.fit.param <- function(processes, states, type.sojourn, distr, init.estim = init.estim, cens.end, cens.beg) {
  
  states <- processes$states
  s <- processes$s
  kmax <- processes$kmax
  counting <- processes$counting
  
  
  if (type.sojourn == "fij") {
    
    if (cens.end) {# We can't decompose the log-likelihood as a sum of optimization problems
      
      ptrans <- matrix(0, nrow = s, ncol = s)
      param <- array(data = NA, dim = c(s, s, 2))
      
      for (i in 1:s) {
        estparam <- .fit.param.fij.endcensoring(counting, i, s, kmax, distr, cens.beg)
        ptrans[i, ] <- estparam$ptrans
        param[i, , ] <- estparam$param
      }
      
    } else {
      
      # Estimation of the transition matrix - component #1
      ptrans <- counting$Nij / matrix(data = counting$Ni, nrow = nrow(counting$Nij), ncol = ncol(counting$Nij))
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
    ptrans <- counting$Nij / matrix(data = counting$Ni, nrow = nrow(counting$Nij), ncol = ncol(counting$Nij))
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
    
    if (cens.end) {# We can't decompose the log-likelihood as a sum of optimization problems
      
      estparam <- .fit.param.fj.endcensoring(counting, s, kmax, distr, cens.beg)
      ptrans <- estparam$ptrans
      param <- estparam$param
      
    } else {
      
      # Estimation of the transition matrix - component #1
      ptrans <- counting$Nij / matrix(data = counting$Ni, nrow = nrow(counting$Nij), ncol = ncol(counting$Nij))
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
    ptrans <- counting$Nij / matrix(data = counting$Ni, nrow = nrow(counting$Nij), ncol = ncol(counting$Nij))
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
  
  # Renormalize ptrans for potential numerical issues
  ptrans <- .normalizePtrans(ptrans)
  
  estimate <-
    smmparametric(
      states = states,
      init = rep.int(x = 1 / s, times = s),
      ptrans = ptrans,
      type.sojourn = type.sojourn,
      distr = distr,
      param = param,
      cens.beg = cens.beg,
      cens.end = cens.end
    )
  
  # Initial distribution
  if (is.vector(init.estim) & length(init.estim) == 1) {
    if (init.estim == "mle") {
      estimate$init <- counting$Nstarti / sum(counting$Nstarti)
    } else if (init.estim == "limit") {
      q <- getKernel(estimate, kmax)
      estimate$init <- .limitDistribution(q = q, ptrans = estimate$ptrans)
    } else if (init.estim == "freq") {
      estimate$init <- counting$Ni / counting$N
    } else if (init.estim == "unif") {
      estimate$init <- rep.int(x = 1 / s, times = s)
    } else {
      stop("'init.estim' must be equal to \"mle\", \"limit\", \"freq\" or \"unif\".
           'init.estim' can also be a vector of length s for custom initial distribution")
    }
  } else {
    if (!(is.numeric(init.estim) & !anyNA(init.estim) & is.vector(init.estim) & length(init.estim) == s)) {
      stop("'init.estim' is not a numeric vector of length s")
    }
    
    if (!(all(init.estim >= 0) & all(init.estim <= 1))) {
      stop("Probabilities in 'init.estim' must be between [0, 1]")
    }
    
    if (!((sum(init.estim) >= 1 - sqrt(.Machine$double.eps)) | (sum(init.estim) <= 1 + sqrt(.Machine$double.eps)))) {
      stop("The sum of 'init.estim' is not equal to one")
    }
    
    estimate$init <- init.estim
  }
  
  estimate$init <- estimate$init / sum(estimate$init)
  names(estimate$init) <- colnames(estimate$ptrans)
  
  return(estimate)
  
}
