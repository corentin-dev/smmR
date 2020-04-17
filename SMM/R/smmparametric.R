#' Parametric semi-Markov model specification
#'
#' @description Creates a parametric model specification for a semi-Markov model.
#'
#' @details This function creates a semi-Markov model object in the parametric 
#' case, taking into account the type of sojourn time and the censoring 
#' described in references. For the parametric specification, several discrete 
#' distributions are considered (see below).
#'
#' The difference between the Markov model and the semi-Markov model concerns 
#' the modelisation of the sojourn time. With a Markov chain, the sojourn time 
#' distribution is modeled by a Geometric distribution (in discrete time). 
#' With a semi-Markov chain, the sojourn time can be any arbitrary distribution.
#' In this package, the available distribution for a semi-Markov model are :
#'  \itemize{
#'    \item Uniform: \eqn{f(x) = 1/n} for \eqn{a \le x \le b}, with \eqn{n = b-a+1};
#'    \item Geometric: \eqn{f(x) = \theta (1-\theta)^x} for \eqn{x = 0, 1, 2,\ldots,n}, \eqn{0 < \theta \le 1}, with \eqn{n > 0} and \eqn{\theta} is the probability of success;
#'    \item Poisson: \eqn{f(x) = (\lambda^x exp(-\lambda))/x!} for \eqn{x = 0, 1, 2,\ldots,n}, with \eqn{n > 0} and \eqn{\lambda > 0};
#'    \item Discrete Weibull of type 1: \eqn{f(x)=q^{(x-1)^{\beta}}-q^{x^{\beta}}, x=1,2,3,\ldots,n}, with \eqn{n > 0}, \eqn{q} is the first parameter and \eqn{\beta} is the second parameter;
#'    \item Negative binomial: \eqn{f(x)=\Gamma(x+\alpha)/(\Gamma(\alpha) x!) (\alpha/(\alpha+\mu))^{\alpha} (\mu/(\alpha+\mu))^x}, for \eqn{x = 0, 1, 2,\ldots,n}, \eqn{n > 0,\ \Gamma } is the Gamma function, \eqn{\alpha} is the parameter of overdispersion and \eqn{\mu} is the mean;
#'    \item Non-parametric.
#'  }
#' We define :
#'  \itemize{
#'    \item the semi-Markov kernel \eqn{q_{ij}(k) = P( J_{m+1} = j, T_{m+1} - T_{m} = k | J_{m} = i )};
#'    \item the transition matrix \eqn{(p_{trans}(i,j))_{i,j} \in E} of the embedded Markov chain \eqn{J = (J_m)_m}, \eqn{p_{trans}(i,j) = P( J_{m+1} = j | J_m = i )};
#'    \item the initial distribution \eqn{\mu_i = P(J_1 = i) = P(Y_1 = i)}, \eqn{i \in 1, 2, \dots, S};
#'    \item the conditional sojourn time distributions \eqn{(f_{ij}(k))_{i,j} \in E,\ k \in N ,\ f_{ij}(k) = P(T_{m+1} - T_m = k | J_m = i, J_{m+1} = j )}, f is specified by the argument "param" in the parametric case and by "laws" in the non-parametric case.
#'  }
#'
#' In this package we can choose differents types of sojourn time. 
#' Four options are available for the sojourn times:
#' \itemize{
#'   \item depending on the present state and on the next state (\eqn{f_{ij}});
#'   \item depending only on the present state (\eqn{f_{i}});
#'   \item depending only on the next state (\eqn{f_{j}});
#'   \item depending neither on the current, nor on the next state (\eqn{f}).
#' }
#' 
#' If  `type.sojourn = "fij"`, `distr` is a matrix of size SxS (e.g., if the 
#' row 1 of the 2nd column is `"pois"`, that is to say we go from the first 
#' state to the second state following a Poisson distribution).
#' If `type.sojourn = "fi"` or `"fj"`, `distr` must be a vector (e.g., if the 
#' first element of vector is `"geom"`, that is to say we go from the first 
#' state to any state according to a Geometric distribution).
#' If `type.sojourn = "f"`, `distr` must be one of `"unif"`, `"geom"`, `"pois"`, 
#' `"dweibull"`, `"nbinom"` (e.g., if `distr` is equal to `"nbinom"`, that is 
#' to say that the sojourn times when going from any state to any state follows 
#' a Negative Binomial distribution).
#' For the non-parametric case, `distr` is equal to `"nonparametric"` whatever 
#' type of sojourn time given.
#' 
#' If the sequence is censored at the beginning and at the end, `cens.beg` 
#' must be equal to `TRUE` and `cens.end` must be equal to `TRUE` too. 
#' All the sequences must be censored in the same way.
#'
#' @param E Vector of state space of length S.
#' @param init Vector of initial distribution of length S.
#' @param ptrans Matrix of transition probabilities of the embedded Markov chain 
#'   \eqn{J=(J_m)_{m}} of size SxS.
#' @param type.sojourn Type of sojourn time (for more explanations, see Details).
#' @param distr
#'   \itemize{
#'     \item Matrix of distributions of size SxS if `type.sojourn = "fij"`;
#'     \item Vector of distributions of size S if `type.sojourn = "fi"` or `"fj`;
#'     \item A distribution if `type.sojourn = "f"`.
#'   }
#'   where the distributions to be used can be one of `unif`, `geom`, `pois`, `dweibull` or `nbinom`.
#' @param param Parameters of sojourn time distributions:
#'   \itemize{
#'     \item Array of distribution parameters of size SxSx2 
#'       (2 corresponds to the maximal number of distribution parameters) if `type.sojourn = "fij"`;
#'     \item Matrix of distribution parameters of size Sx2 if `type.sojourn = "fi"` or `"fj"`;
#'     \item Vector of distribution parameters of length 2 if `type.sojourn = "f"`.
#'   }
#'   
#'  When parameters/values are not necessary (e.g. the Poisson distribution has 
#'  only one parameter that is \eqn{\lambda}, leave the value `NA` for the 
#'  second parameter in the argument `param`).
#' @param cens.beg Optional. A logical value indicating whether or not 
#'   sequences are censored at the beginning.
#' @param cens.end Optional. A logical value indicating whether or not 
#'   sequences are censored at the end.
#' @return Returns an object of class [smmparametric][smmparametric].
#' 
#' 
#' @seealso [simulate], [fitsemimarkovmodel], [smmnonparametric]
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
#' smm1
#' 
smmparametric <- function(E, init, ptrans, type.sojourn = c("fij", "fi", "fj", "f"), 
                          distr, param, cens.beg = FALSE, cens.end = FALSE) {
  
  #############################
  # Checking parameter E
  #############################
  
  S <- length(E)
  
  if (!(is.vector(E) && (length(unique(E)) == S))) {
    stop("The state space E is not a vector of unique elements")
  }
  
  #############################
  # Checking parameter init
  #############################
  
  if (!(is.vector(init) && (length(init) == S))) {
    stop("init is not a vector of length S")
  }
  
  if (!(all(init >= 0) && all(init <= 1))) {
    stop("Probabilities in init must be between [0, 1]")
  }
  
  if (!(sum(init) == 1)) {
    stop("The sum of init is not equal to one")
  }
  
  #############################
  # Checking parameter ptrans
  #############################
  
  if (!is.matrix(ptrans)) {
    stop("ptrans is not a matrix")
  }
  
  if (!(all(ptrans >= 0) && all(ptrans <= 1))) {
    stop("Probabilities in ptrans must be between [0, 1]")
  }
  
  if (!((dim(ptrans)[1] == S) && (dim(ptrans)[2] == S))) {
    stop("The size of the matrix ptrans must be equal to SxS")
  }
  
  if (!all(apply(ptrans, 1, sum) == 1)) {
    stop("ptrans is not a stochastic matrix (column sums accross rows must be equal to one for each row)")
  }
  
  if (!all(diag(ptrans) == 0)) {
    stop("All the diagonal elements of ptrans must be equal to 0 since transitions to the same state are not allowed")
  }
  
  #############################
  # Checking parameter type.sojourn
  #############################
  
  type.sojourn <- match.arg(type.sojourn)
  
  #############################
  # Checking parameters distr and param
  #############################
  
  if (type.sojourn == "fij" && !(is.matrix(distr) && (is.array(param) && !is.matrix(param)))) {
    stop("distr must be a matrix of size SxS and param must be an array of size SxSx2 since type.sojourn == \"fij\"")
  }
  
  if ((type.sojourn == "fi" | type.sojourn == "fj") && !(is.vector(distr) && is.matrix(param))) {
    stop("distr must be a vector of length S and param must be a matrix of size Sx2 since type.sojourn == \"fi\" or \"fj\"")
  }
  
  if (type.sojourn == "f" && !((length(distr) == 1) && is.vector(param))) {
    stop("distr must be one element and param must be a vector of length 2 since type.sojourn == \"f\"")
  }
  
  
  if (type.sojourn == "fij" && !((dim(distr)[1] == S && dim(distr)[2] == S) && (dim(param)[1] == S && dim(param)[2] == S))) {
    stop("distr must be a matrix of size SxS and param must be an array of size SxSx2 since type.sojourn == \"fij\"")
  }
  
  if ((type.sojourn == "fi" | type.sojourn == "fj") && !((length(distr) == S) && (dim(param)[1] == S && dim(param)[2] == 2))) {
    stop("distr must be a vector of length S and param must be a matrix of size Sx2 since type.sojourn == \"fi\" or \"fj\"")
  }
  
  
  distrib.vec <- c("unif", "geom", "pois", "dweibull", "nbinom", NA)
  if (!all(distr %in% distrib.vec)) {
    stop("The specified distributions must be either ", paste(distrib.vec, collapse = ", "), 
         ".\n Incorrect distribution(s) found in distr: ", paste(as.character(distr)[!(distr %in% distrib.vec)], collapse = ", "))
  }
  
  if (!all(param >= 0 | is.na(param))) {
    stop("Every element of param must be positive")
  }
  
  if (type.sojourn == "fij") {
    if (!(all(is.na(diag(distr))))) {
      stop("All the diagonal elements of distr must be equal to NA since transitions to the same state are not allowed")
    }
    
    if (!(all(is.na(diag(param[, , 1]))) && all(is.na(diag(param[, , 2]))))) {
      stop("All the diagonal elements of param must be equal to NA since transitions to the same state are not allowed")
    }
    
    if (!all(!is.na(distr[row(distr) != col(distr)]))) {
      stop("All the non-diagonal elements of distr must be specified. Found NAs values.")
    }
  }
  
  #############################
  # Checking parameter cens.beg and cens.end
  #############################
  
  if (!(is.logical(cens.beg) && is.logical(cens.end))) {
    stop("cens.beg and cens.end must be TRUE or FALSE")
  }
  
  
  ans <-
    list(
      E = E,
      S = S,
      init = init,
      type.sojourn = type.sojourn,
      ptrans = ptrans,
      distr = distr,
      param = param,
      cens.beg = cens.beg,
      cens.end = cens.end
    )
  
  class(ans) <- c("smm", "smmparametric")
  
  return(ans)
}

is.smmparametric <- function(x) {
  inherits(x, "smmparametric")
}

.get.f.smmparametric <- function(x, Kmax) {
  
  S <- x$S
  
  if (x$type.sojourn == "fij") {
    param1 <- x$param[, , 1]
    param2 <- x$param[, , 2]
    f <- matrix(0, nrow = S * S, ncol = Kmax)
  } else if (x$type.sojourn == "fj") {
    param1 <- x$param[, 1]
    param2 <- x$param[, 2]
    f <- matrix(0, nrow = S, ncol = Kmax)
  } else if (x$type.sojourn == "fi") {
    param1 <- x$param[, 1]
    param2 <- x$param[, 2]
    f <- matrix(0, nrow = S, ncol = Kmax)
  } else {
    param1 <- x$param[1]
    param2 <- x$param[2]
    f <- matrix(0, nrow = 1, ncol = Kmax)
  }
  
  if ("dweibull" %in% x$distr) {
    indices <- which(x$distr == "dweibull")
    for (j in indices) {
      f[j, ] <- ddweibull(1:Kmax, q = param1[j], beta = param2[j], zero = FALSE)
    }
  }
  if ("geom" %in% x$distr) {
    indices <- which(x$distr == "geom")
    for (j in indices) {
      f[j, ] <- dgeom(0:(Kmax - 1), prob = param1[j])
    }
  }
  if ("nbinom" %in% x$distr) {
    indices <- which(x$distr == "nbinom")
    for (j in indices) {
      f[j, ] <- dnbinom(0:(Kmax - 1), size = param1[j], mu = param2[j])
    }
  }
  if ("pois" %in% x$distr) {
    indices <- which(x$distr == "pois")
    for (j in indices) {
      f[j, ] <- dpois(0:(Kmax - 1), lambda = param1[j])
    }
  }
  if ("unif" %in% x$distr) {
    indices <- which(x$distr == "unif")
    for (j in indices) {
      f[j, ] <- sapply(1:Kmax, function(k) ifelse(k <= x$param[j], 1 / x$param[j], 0))
    }
  }
  
  if (x$type.sojourn == "fij") {
    fijk <- array(f, c(S, S, Kmax))
  } else if (x$type.sojourn == "fi") {
    f <- rep(as.vector(t(f)), each = S)
    fmat <- matrix(f, nrow = Kmax, ncol = S * S, byrow = TRUE)
    fk <- array(as.vector(t(fmat)), c(S, S, Kmax))
    fijk <- apply(X = fk, MARGIN =  c(1, 3), FUN =  t)
  } else if (x$type.sojourn == "fj") {
    f <- rep(as.vector(t(f)), each = S)
    fmat <- matrix(f, nrow = Kmax, ncol = S * S, byrow = TRUE)
    fijk <- array(as.vector(t(fmat)), c(S, S, Kmax))
  } else {
    f <- rep(f, each = S * S)
    fmat <- matrix(f, nrow = Kmax, ncol = S * S, byrow = TRUE)
    fijk <- array(as.vector(t(fmat)), c(S, S, Kmax))
  }
  
  return(fijk)
  
}

.get.q.smmparametric <- function(x, Kmax) {
  
  S <- x$S
  
  fijk <- .get.f(x, Kmax)
  q <- array(x$ptrans, c(S, S, Kmax)) * fijk
  
  return(q)
  
}

.getKpar.smmparametric <- function(x) {
  
  distr <- x$distr
  
  nbDweibull <- length(which(distr == "dweibull"))
  nbGeom <- length(which(distr == "geom"))
  nbNbinom <- length(which(distr == "nbinom"))
  nbPois <- length(which(distr == "pois"))
  nbUnif <- length(which(distr == "unif"))
  
  Kpar <- 2 * nbDweibull + nbGeom + 2 * nbNbinom + nbPois + nbUnif
  
  return(Kpar)
}