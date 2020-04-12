#' Simulates k-th order Markov chains
#'
#' @description Simulates k-th order Markov chains.
#'
#' @details If `nsim` is a single integer then a chain of that length is 
#'   produced. If `nsim` is a vector of integers, then `length(nsim)` 
#'   sequences are generated with respective lengths.
#'
#' @param object An object of class [markovmodel].
#' @param nsim An integer or vector of integers (for multiple sequences) 
#'   specifying the length of the sequence(s).
#' @param seed `seed` for the random number generator.
#' @param ... further arguments passed to or from other methods.
#' @return A list of vectors representing the sequences.
#' 
#' 
#' @seealso [fitmarkovmodel], [smmnonparametric], [smmparametric], [fitsemimarkovmodel]
#' @export
#'
simulate.markovmodel <- function(object, nsim = 1, seed = NULL, ...) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  S <- length(object$E)
  out <- list()
  nbseq <- length(nsim)
  
  for (n in 1:nbseq) {
    y <- rep.int(NA, nsim[n])
    
    y[1:object$k] <- sample(x = object$E, size = object$k, replace = TRUE, prob = object$init)
    
    for (i in 1:(nsim[n] - object$k)) {
      ind <- which(object$E == y[i + object$k - 1])
      if (object$k > 1) {
        for (j in (object$k - 2):0) {
          ind <- ind + S ^ (j + 1) * (which(object$E == y[i + j]) - 1)
        }
      }
      y[i + object$k] <- sample(object$E, 1, prob = object$ptrans[ind, ])
    }
    out[[n]] <- y
    
  }
  
  return(out)
  
}
