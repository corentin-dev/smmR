#' Simulates semi-Markov chains
#'
#' @description Simulates sequences from a parametric semi-Markov model.
#'
#' @details If `nsim` is a single integer then a chain of that length is 
#'   produced. If `nsim` is a vector of integers, then `length(nsim)` 
#'   sequences are generated with respective lengths.
#'
#' @param object An object of class [smmparametric].
#' @param nsim An integer or vector of integers (for multiple sequences) 
#'   specifying the length of the sequence(s).
#' @param seed `seed` for the random number generator.
#' @param ... further arguments passed to or from other methods.
#' @return A list of vectors representing the sequences.
#' 
#' 
#' @seealso [smmnonparametric], [fitsemimarkovmodel]
#' @export
#'
#' @examples 
#' E <- c("a", "c", "g", "t")
#' S <- length(E)
#' 
#' # Creation of the initial distribution
#' vect.init <- c(1 / 4, 1 / 4, 1 / 4, 1 / 4)
#' 
#' # Creation of transition matrix
#' pij <- matrix(c(0, 0.2, 0.5, 0.3, 
#'                 0.2, 0, 0.3, 0.5, 
#'                 0.3, 0.5, 0, 0.2, 
#'                 0.4, 0.2, 0.4, 0), 
#'               ncol = S, byrow = TRUE)
#' 
#' # Creation of the distribution matrix
#' 
#' distr.matrix <- matrix(c(NA, "pois", "geom", "nbinom", 
#'                          "geom", NA, "pois", "dweibull",
#'                          "pois", "pois", NA, "geom", 
#'                          "pois", "geom", "geom", NA), 
#'                        nrow = S, ncol = S, byrow = TRUE)
#' 
#' # Creation of an array containing the parameters
#' param1.matrix <- matrix(c(NA, 2, 0.4, 4, 
#'                           0.7, NA, 5, 0.6, 
#'                           2, 3, NA, 0.6, 
#'                           4, 0.3, 0.4, NA), 
#'                         nrow = S, ncol = S, byrow = TRUE)
#' 
#' param2.matrix <- matrix(c(NA, NA, NA, 2, 
#'                           NA, NA, NA, 0.8, 
#'                           NA, NA, NA, NA, 
#'                           NA, NA, NA, NA), 
#'                         nrow = S, ncol = S, byrow = TRUE)
#' 
#' param.array <- array(c(param1.matrix, param2.matrix), c(S, S, 2))
#' 
#' # Specify the semi-Markov model
#' smm1 <- smmparametric(E = E, init = vect.init, ptrans = pij, 
#'                       type.sojourn = "fij", distr = distr.matrix, 
#'                       param = param.array)
#' 
#' seq1 <- simulate(object = smm1, nsim = c(1000, 10000, 2000), seed = 100)
#' seq1[[1]][1:15]
#' 
simulate.smmparametric <- function(object, nsim = 1, seed = NULL, ...) {
  
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
        
        if (object$distr[which(J[i] == object$E), which(J[i + 1] == object$E)] == "unif") {
          
          Kmax <- object$param[which(J[i] == object$E), which(J[i + 1] == object$E), 1]
          k <- sample(1:Kmax, 1)
          
        } else if (object$distr[which(J[i] == object$E), which(J[i + 1] == object$E)] == "geom") {
          
          k <- rgeom(1, object$param[which(J[i] == object$E), which(J[i + 1] == object$E), 1]) + 1
          
        } else if (object$distr[which(J[i] == object$E), which(J[i + 1] == object$E)] == "pois") {
          
          k <- rpois(1, object$param[which(J[i] == object$E), which(J[i + 1] == object$E), 1]) + 1 # no sojourn time equal to zero so we shift the Poisson distribution
          
        } else if (object$distr[which(J[i] == object$E), which(J[i + 1] == object$E)] == "dweibull") {
          
          k <- rdweibull(1, object$param[which(J[i] == object$E), which(J[i + 1] == object$E), 1], object$param[which(J[i] == object$E), which(J[i + 1] == object$E), 2], zero = FALSE)
          
        } else if (object$distr[which(J[i] == object$E), which(J[i + 1] == object$E)] == "nbinom") {
          
          k <- rnbinom(n = 1, size = object$param[which(J[i] == object$E), which(J[i + 1] == object$E), 1], mu = object$param[which(J[i] == object$E), which(J[i + 1] == object$E), 2]) + 1
          
        }
        
      } else if (object$type.sojourn == "fi") {
        
        if (object$distr[which(J[i] == object$E)] == "unif") {
          
          Kmax <- object$param[which(J[i] == object$E), 1]
          k <- sample(1:Kmax, 1)
          
        } else if (object$distr[which(J[i] == object$E)] == "geom") {
          
          k <- rgeom(1, object$param[which(J[i] == object$E), 1]) + 1
          
        } else if (object$distr[which(J[i] == object$E)] == "pois") {
          
          k <- rpois(1, object$param[which(J[i] == object$E), 1]) + 1
          
        } else if (object$distr[which(J[i] == object$E)] == "dweibull") {
          
          k <- rdweibull(1, object$param[which(J[i] == object$E), 1], object$param[which(J[i] == object$E), 2], zero = FALSE)
          
        } else if (object$distr[which(J[i] == object$E)] == "nbinom") {
          
          k <- rnbinom(n = 1, size = object$param[which(J[i] == object$E), 1], mu = object$param[which(J[i] == object$E), 2]) + 1
          
        }
        
      } else if (object$type.sojourn == "fj") {
        
        if (object$distr[which(J[i + 1] == object$E)] == "unif") {
          
          Kmax <- object$param[which(J[i + 1] == object$E), 1]
          k <- sample(1:Kmax, 1)
          
        } else if (object$distr[which(J[i + 1] == object$E)] == "geom") {
          
          k <- rgeom(1, object$param[which(J[i + 1] == object$E), 1]) + 1
          
        } else if (object$distr[which(J[i + 1] == object$E)] == "pois") {
          
          k <- rpois(1, object$param[which(J[i + 1] == object$E), 1]) + 1
          
        } else if (object$distr[which(J[i + 1] == object$E)] == "dweibull") {
          
          k <- rdweibull(1, object$param[which(J[i + 1] == object$E), 1], object$param[which(J[i + 1] == object$E), 2], zero = FALSE)
          
        } else if (object$distr[which(J[i + 1] == object$E)] == "nbinom") {
          
          k <- rnbinom(n = 1, size = object$param[which(J[i + 1] == object$E), 1], mu = object$param[which(J[i + 1] == object$E), 2]) + 1
          
        }
        
      } else {# f case
        
        if (object$distr == "unif") {
          
          Kmax <- object$param[1]
          k <- sample(1:Kmax, 1)
          
        } else if (object$distr == "geom") {
          
          k <- rgeom(1, object$param[1]) + 1 # we shift
          
        } else if (object$distr == "pois") {
          
          k <- rpois(1, object$param[1]) + 1 # we shift
          
        } else if (object$distr == "dweibull") {
          
          k <- rdweibull(1, object$param[1], object$param[2], zero = FALSE)
          
        } else if (object$distr == "nbinom") {
          
          k <- rnbinom(n = 1, size = object$param[1], mu = object$param[2]) + 1 # we shift
          
        }
        
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
    
    out[[m]] <- y
    
    
  }
  
  return(out)
  
}