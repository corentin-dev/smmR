#' Simulates semi-Markov chains
#'
#' @description Simulates sequences from a non-parametric semi-Markov model.
#'
#' @details If `nsim` is a single integer then a chain of that length is 
#'   produced. If `nsim` is a vector of integers, then `length(nsim)` 
#'   sequences are generated with respective lengths.
#'
#' @param object An object of class [smmnonparametric].
#' @param nsim An integer or vector of integers (for multiple sequences) 
#'   specifying the length of the sequence(s).
#' @param seed `seed` for the random number generator.
#' @param ... further arguments passed to or from other methods.
#' @return A list of vectors representing the sequences.
#' 
#' 
#' @seealso [smmparametric], [fitsemimarkovmodel]
#' @export
#'
#' @examples
#' states <- c("a", "c", "g", "t")
#' s <- length(states)
#' 
#' # Creation of the initial distribution
#' vect.init <- c(1 / 4, 1 / 4, 1 / 4, 1 / 4)
#' 
#' # Creation of transition matrix
#' pij <- matrix(c(0, 0.2, 0.5, 0.3, 
#'                 0.2, 0, 0.3, 0.5, 
#'                 0.3, 0.5, 0, 0.2, 
#'                 0.4, 0.2, 0.4, 0), 
#'               ncol = s, byrow = TRUE)
#' 
#' # Creation of a matrix corresponding to the 
#' # conditional sojourn time distributions
#' Kmax <- 6
#' nparam.matrix <- matrix(c(0.2, 0.1, 0.3, 0.2, 
#'                           0.2, 0, 0.4, 0.2, 
#'                           0.1, 0, 0.2, 0.1, 
#'                           0.5, 0.3, 0.15, 0.05, 
#'                           0, 0, 0.3, 0.2, 
#'                           0.1, 0.2, 0.2, 0), 
#'                         nrow = s, ncol = Kmax, byrow = TRUE)
#' 
#' smm2 <- smmnonparametric(states = states, init = vect.init, ptrans = pij, 
#'                          type.sojourn = "fj", distr = nparam.matrix)
#'
#' seq2 <- simulate(object = smm2, nsim = c(1000, 10000, 2000), seed = 100)
#' seq2[[1]][1:15]
#'
simulate.smmnonparametric <- function(object, nsim = 1, seed = NULL, ...) {
  
  ###########################################################
  ###########################################################
  # The algorithm used to simulate the sequences is the following:
  # 
  # 1. Set k = 0, S_{0} = 0 and sample J_{0} from the initial distribution \alpha;
  # 2. Sample the random variable J \sim p(J_{k} , .) and set J_{k+1} = J(\omega);
  # 3. Sample the random variable X \sim F_{J_{k} J_{k+1}}(.)
  # 4. Set S_{k+1} = S_{k} + X;
  # 5. If S_{k+1} >= M, then end;
  # 6. Else, set k= k + 1 and continue to step 2.
  # 
  ###########################################################
  ###########################################################
  
  if (!is.null(seed)) {
    set.seed(seed)  
  }
  
  sequences <- list()
  nbseq <- length(nsim)
  
  for (m in 1:nbseq) {
    
    J <- c()
    T <- c()
    J[1] <- sample(object$states, 1, prob = object$init)
    
    i <- 1
    t <- 1
    
    while (t <= nsim[m]) {
      
      J[i + 1] <- sample(object$states, 1, prob = object$ptrans[which(object$states == J[i]), ])
      
      if (object$type.sojourn == "fij") {
        
        Kmax <- dim(object$distr)[3]
        k <- sample(1:Kmax, 1, prob = object$distr[which(J[i] == object$states), which(J[i + 1] == object$states), ])
        
      } else if (object$type.sojourn == "fi") {
        
        Kmax <- dim(object$distr)[2]
        k <- sample(1:Kmax, 1, prob = object$distr[which(J[i] == object$states), ])
        
      } else if (object$type.sojourn == "fj") {
        
        Kmax <- dim(object$distr)[2]
        k <- sample(1:Kmax, 1, prob = object$distr[which(J[i + 1] == object$states), ])
        
      } else {
        
        Kmax <- length(object$distr)
        k <- sample(1:Kmax, 1, prob = object$distr)
        
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
      
      y <- .getSeq(J, T)
      y <- y[Nl:(t - 1 - Nl)]
      
    } else if (object$cens.beg == FALSE && object$cens.end == TRUE) {# First time is a Jump Time
      
      y <- .getSeq(J, T)
      y <- y[1:nsim[m]]
      
    } else if (object$cens.beg == 1 && object$cens.end == 0) {
      
      l <- t - nsim[m]
      y <- .getSeq(J, T)
      y <- y[l:(t - 1)]
      
    } else {# First and last times are jump times
      
      y <- .getSeq(J, T)
      
    }
    
    sequences[[m]] <- y
  }
  
  
  return(sequences)
  
}
