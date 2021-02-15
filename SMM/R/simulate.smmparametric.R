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
#' @seealso [smmparametric], [smmnonparametric], [fitsemimarkovmodel]
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
#' # Creation of the distribution matrix
#' 
#' distr.matrix <- matrix(c(NA, "pois", "geom", "nbinom", 
#'                          "geom", NA, "pois", "dweibull",
#'                          "pois", "pois", NA, "geom", 
#'                          "pois", "geom", "geom", NA), 
#'                        nrow = s, ncol = s, byrow = TRUE)
#' 
#' # Creation of an array containing the parameters
#' param1.matrix <- matrix(c(NA, 2, 0.4, 4, 
#'                           0.7, NA, 5, 0.6, 
#'                           2, 3, NA, 0.6, 
#'                           4, 0.3, 0.4, NA), 
#'                         nrow = s, ncol = s, byrow = TRUE)
#' 
#' param2.matrix <- matrix(c(NA, NA, NA, 0.6, 
#'                           NA, NA, NA, 0.8, 
#'                           NA, NA, NA, NA, 
#'                           NA, NA, NA, NA), 
#'                         nrow = s, ncol = s, byrow = TRUE)
#' 
#' param.array <- array(c(param1.matrix, param2.matrix), c(s, s, 2))
#' 
#' # Specify the semi-Markov model
#' smm1 <- smmparametric(states = states, init = vect.init, ptrans = pij, 
#'                       type.sojourn = "fij", distr = distr.matrix, 
#'                       param = param.array)
#' 
#' seq1 <- simulate(object = smm1, nsim = c(1000, 10000, 2000), seed = 100)
#' seq1[[1]][1:15]
#' 
simulate.smmparametric <- function(object, nsim = 1, seed = NULL, ...) {
  
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
  
  # Preparation of the parameters and the distributions matrix to ease the sampling process
  param1 <- rep.int(x = NA, times = object$s)
  param2 <- rep.int(x = NA, times = object$s)
  distributions <- array(data = NA, dim = c(object$s, object$s))
  if (object$type.sojourn == "fij") {
    param1 <- object$param[, , 1]
    param2 <- object$param[, , 2]
    distributions <- object$distr
  } else if (object$type.sojourn == "fj") {
    param1 <- t(matrix(data = object$param[, 1], nrow = object$s, ncol = object$s))
    param2 <- t(matrix(data = object$param[, 2], nrow = object$s, ncol = object$s))
    distributions <- t(matrix(data = object$distr, nrow = object$s, ncol = object$s))
  } else if (object$type.sojourn == "fi") {
    param1 <- matrix(data = object$param[, 1], nrow = object$s, ncol = object$s)
    param2 <- matrix(data = object$param[, 2], nrow = object$s, ncol = object$s)
    distributions <- matrix(data = object$distr, nrow = object$s, ncol = object$s)
  } else {
    param1 <- matrix(data = object$param[1], nrow = object$s, ncol = object$s)
    param2 <- matrix(data = object$param[2], nrow = object$s, ncol = object$s)
    distributions <- matrix(data = object$distr, nrow = object$s, ncol = object$s)
  }
  
  distrib <- matrix(data = NA, nrow = nrow(distributions), ncol = ncol(distributions))
  
  distrib[distributions == "unif"] <- 0
  distrib[distributions == "geom"] <- 1
  distrib[distributions == "pois"] <- 2
  distrib[distributions == "dweibull"] <- 3
  distrib[distributions == "nbinom"] <- 4
  
  sequences <- simulateParam(seed, nsim, object$init, object$ptrans, distrib, param1,
                           param2, censBeg = object$cens.beg, censEnd = object$cens.end)
  
  sequences <- lapply(sequences, function(x) object$states[x])
  
  return(sequences)
  
}
