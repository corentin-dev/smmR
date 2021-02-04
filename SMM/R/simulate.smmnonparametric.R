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
#' @seealso [smmnonparametric], [smmparametric], [fitsemimarkovmodel]
#' 
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
  
  if (is.null(seed)) {
    seed <- as.numeric(Sys.time())
  }
  
  # Preparation of distribution matrix to ease the sampling process
  distribution <- array(data = NA, dim = c(object$s, object$s, object$kmax))
  if (object$type.sojourn == "fij") {
    distribution <- object$distr
  } else if (object$type.sojourn == "fj") {
    distribution <- aperm(a = array(data = object$distr, dim = c(object$s, object$kmax, object$s)), perm = c(3, 1, 2))
  } else if (object$type.sojourn == "fi") {
    distribution <- aperm(a = array(data = object$distr, dim = c(object$s, object$kmax, object$s)), perm = c(1, 3, 2))
  } else {
    distribution <- aperm(a = array(data = object$distr, dim = c(object$kmax, object$s, object$s)), perm = c(2, 3, 1))
  }
  
  sequences <- simulateNonParam(seed, nsim, object$init, object$ptrans, distribution, 
                                censBeg = object$cens.beg, censEnd = object$cens.end)
  
  sequences <- lapply(sequences, function(x) object$states[x])
  
  return(sequences)
  
}
