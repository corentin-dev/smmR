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
#' @seealso [markovmodel], [fitmarkovmodel]
#' 
#' @export
#'
#' @examples 
#' states <- c("a", "c", "g", "t")
#' s <- length(states)
#' k <- 2
#' vect.init <- rep.int(1 / s ^ k, s ^ k)
#' p <- matrix(0.25, nrow = s ^ k, ncol = s)
#' 
#' # Specify the Markov model
#' markov1 <- markovmodel(states = states, init = vect.init, ptrans = p, k = k)
#' 
#' seq1 <- simulate(object = markov1, nsim = c(1000, 10000, 2000), seed = 150)
#' seq1[[1]][1:15]
#' 
simulate.markovmodel <- function(object, nsim = 1, seed = NULL, ...) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  s <- length(object$states)
  out <- list()
  nbseq <- length(nsim)
  
  for (n in 1:nbseq) {
    y <- rep.int(NA, nsim[n])
    
    # Initial state(s)
    y[1:object$k] <- s2c(sample(x = names(object$init), size = 1, prob = object$init))
    
    for (i in 1:(nsim[n] - object$k)) {
      ind <- which(object$states == y[i + object$k - 1])
      if (object$k > 1) {
        for (j in (object$k - 2):0) {
          ind <- ind + s ^ (j + 1) * (which(object$states == y[i + j]) - 1)
        }
      }
      y[i + object$k] <- sample(object$states, 1, prob = object$ptrans[ind, ])
    }
    out[[n]] <- y
    
  }
  
  return(out)
  
}
