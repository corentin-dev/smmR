simulate.smmnonparametric <- function(object, nsim = 1, seed = NULL, ...) {
  
  # If nsim is a single integer then a SMM of that length is produced. 
  #   If nsim is a vector of integers, then length(nsim) sequences are 
  #   generated with respective lengths.
  
  if (!is.null(seed)) {
    set.seed(seed)  
  }
  
  out <- list()
  nbseq <- length(nsim)
  
  for (m in 1:nbseq) {
    
    J <- c()
    T <- c()
    J[1] <- sample(object$E, 1, prob = object$init)
    
    i <- 1
    t <- 1
    
    while (t <= nsim[m]) {
      
      J[i + 1] <- sample(object$E, 1, prob = object$ptrans[which(object$E == J[i]), ])
      
      if (object$type.sojourn == "fij") {
        
        Kmax <- dim(object$laws)[3]
        k <- sample(1:Kmax, 1, prob = object$laws[which(J[i] == object$E), which(J[i + 1] == object$E), ])
        
      } else if (object$type.sojourn == "fi") {
        
        Kmax <- dim(object$laws)[2]
        k <- sample(1:Kmax, 1, prob = object$laws[which(J[i] == object$E), ])
        
      } else if (object$type.sojourn == "fj") {
        
        Kmax <- dim(object$laws)[2]
        k <- sample(1:Kmax, 1, prob = object$laws[which(J[i + 1] == object$E), ])
        
      } else {
        
        Kmax <- length(object$laws)
        k <- sample(1:Kmax, 1, prob = object$laws)
        
      }
      
      T[i] <- t + k
      t <- T[i]
      i <- i + 1
      
    }
    
    #############################
    # Censoring sequences
    #############################
    if (object$cens.beg == TRUE && object$cens.end == TRUE) {
      
      l <- t - nsim[m]
      n <- nsim[m]
      Nl <- floor(l / 2)
      
      y <- .get.seq(J, T)
      y <- y[Nl:(t - 1 - Nl)]
      
    } else if (object$cens.beg == FALSE && object$cens.end == TRUE) {# First time is a Jump Time
      
      y <- .get.seq(J, T)
      y <- y[1:nsim[m]]
      
    } else if (object$cens.beg == 1 && object$cens.end == 0) {
      
      l <- t - nsim[m]
      y <- .get.seq(J, T)
      y <- y[l:(t - 1)]
      
    } else {# First and last times are jump times
      
      y <- .get.seq(J, T)
      
    }
    
    out[[m]] <- y
  }
  
  
  return(out)
  
}
