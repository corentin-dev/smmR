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
#' the modeling of the sojourn time. With a Markov chain, the sojourn time 
#' distribution is modeled by a Geometric distribution (in discrete time). 
#' With a semi-Markov chain, the sojourn time can be any arbitrary distribution.
#' In this package, the available distribution for a semi-Markov model are :
#'  \itemize{
#'    \item Uniform: \eqn{f(x) = 1/n} for \eqn{a \le x \le b}, with \eqn{n = b-a+1};
#'    \item Geometric: \eqn{f(x) = \theta (1-\theta)^x} for \eqn{x = 0, 1, 2,\ldots,n}, \eqn{0 < \theta < 1}, with \eqn{n > 0} and \eqn{\theta} is the probability of success;
#'    \item Poisson: \eqn{f(x) = (\lambda^x exp(-\lambda))/x!} for \eqn{x = 0, 1, 2,\ldots,n}, with \eqn{n > 0} and \eqn{\lambda > 0};
#'    \item Discrete Weibull of type 1: \eqn{f(x)=q^{(x-1)^{\beta}}-q^{x^{\beta}}, x=1,2,3,\ldots,n}, with \eqn{n > 0}, \eqn{q} is the first parameter and \eqn{\beta} is the second parameter;
#'    \item Negative binomial: \eqn{f(x)=\frac{\Gamma(x+\alpha)}{\Gamma(\alpha) x!} p^{\alpha} (1 - p)^x}, 
#'      for \eqn{x = 0, 1, 2,\ldots,n}, \eqn{\Gamma} is the Gamma function, 
#'      \eqn{\alpha} is the parameter of overdispersion and \eqn{p} is the 
#'      probability of success, \eqn{0 < p < 1};
#'    \item Non-parametric.
#'  }
#' We define :
#'  \itemize{
#'    \item the semi-Markov kernel \eqn{q_{ij}(k) = P( J_{m+1} = j, T_{m+1} - T_{m} = k | J_{m} = i )};
#'    \item the transition matrix \eqn{(p_{trans}(i,j))_{i,j} \in states} of the embedded Markov chain \eqn{J = (J_m)_m}, \eqn{p_{trans}(i,j) = P( J_{m+1} = j | J_m = i )};
#'    \item the initial distribution \eqn{\mu_i = P(J_1 = i) = P(Z_1 = i)}, \eqn{i \in 1, 2, \dots, s};
#'    \item the conditional sojourn time distributions \eqn{(f_{ij}(k))_{i,j} \in states,\ k \in N ,\ f_{ij}(k) = P(T_{m+1} - T_m = k | J_m = i, J_{m+1} = j )},
#'      \eqn{f} is specified by the argument `param` in the parametric case.
#'  }
#'  
#' In this package we can choose different types of sojourn time. 
#' Four options are available for the sojourn times:
#' \itemize{
#'   \item depending on the present state and on the next state (\eqn{f_{ij}});
#'   \item depending only on the present state (\eqn{f_{i}});
#'   \item depending only on the next state (\eqn{f_{j}});
#'   \item depending neither on the current, nor on the next state (\eqn{f}).
#' }
#' 
#' If  `type.sojourn = "fij"`, `distr` is a matrix of dimension \eqn{(s, s)}
#' (e.g., if the row 1 of the 2nd column is `"pois"`, that is to say we go from
#' the first state to the second state following a Poisson distribution).
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
#' If the sequence is censored at the beginning and/or at the end, `cens.beg` 
#' must be equal to `TRUE` and/or `cens.end` must be equal to `TRUE`.
#' All the sequences must be censored in the same way.
#' 
#' @param states Vector of state space of length \eqn{s}.
#' @param init Vector of initial distribution of length \eqn{s}.
#' @param ptrans Matrix of transition probabilities of the embedded Markov 
#'   chain \eqn{J=(J_m)_{m}} of dimension \eqn{(s, s)}.
#' @param type.sojourn Type of sojourn time (for more explanations, see Details).
#' @param distr
#'   \itemize{
#'     \item Matrix of distributions of dimension \eqn{(s, s)} if `type.sojourn = "fij"`;
#'     \item Vector of distributions of length \eqn{s} if `type.sojourn = "fi"` or `"fj`;
#'     \item A distribution if `type.sojourn = "f"`.
#'   }
#'   where the distributions to be used can be one of `unif`, `geom`, `pois`, `dweibull` or `nbinom`.
#' @param param Parameters of sojourn time distributions:
#'   \itemize{
#'     \item Array of distribution parameters of dimension \eqn{(s, s, 2)}
#'       (2 corresponds to the maximal number of distribution parameters) if `type.sojourn = "fij"`;
#'     \item Matrix of distribution parameters of dimension \eqn{(s, 2)} if `type.sojourn = "fi"` or `"fj"`;
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
#' @return Returns an object of class [smmparametric].
#' 
#' @seealso [simulate], [fitsmm], [smmnonparametric]
#' 
#' @references
#' V. S. Barbu, N. Limnios. (2008). Semi-Markov Chains and Hidden Semi-Markov 
#' Models Toward Applications - Their Use in Reliability and DNA Analysis. 
#' New York: Lecture Notes in Statistics, vol. 191, Springer.
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
#' smm1
#' 
smmparametric <- function(states, init, ptrans, type.sojourn = c("fij", "fi", "fj", "f"), 
                          distr, param, cens.beg = FALSE, cens.end = FALSE) {
  
  #############################
  # Checking parameter states
  #############################
  
  s <- length(states)
  
  if (!(is.vector(states) & (length(unique(states)) == s))) {
    stop("The state space 'states' is not a vector of unique elements")
  }
  
  #############################
  # Checking parameter init
  #############################
  
  if (!(is.numeric(init) & !anyNA(init) & is.vector(init) & length(init) == s)) {
    stop("'init' is not a numeric vector of length s")
  }
  
  if (!(all(init >= 0) & all(init <= 1))) {
    stop("Probabilities in 'init' must be between [0, 1]")
  }
  
  if (!((sum(init) >= 1 - .Machine$double.eps) | (sum(init) <= 1 + .Machine$double.eps))) {
    stop("The sum of 'init' is not equal to one")
  }
  
  #############################
  # Checking parameter ptrans
  #############################
  
  if (!(is.numeric(ptrans) & !anyNA(ptrans) & is.matrix(ptrans))) {
    stop("'ptrans' is not a matrix with numeric values")
  }
  
  if (!((dim(ptrans)[1] == s) & (dim(ptrans)[2] == s))) {
    stop("The dimension of the matrix 'ptrans' must be equal to (s, s)")
  }
  
  if (!(all(ptrans >= 0) & all(ptrans <= 1))) {
    stop("Probabilities in 'ptrans' must be between [0, 1]")
  }
  
  if (!all(diag(ptrans) == 0)) {
    stop("All the diagonal elements of 'ptrans' must be equal to 0 since transitions to the same state are not allowed")
  }
  
  if (!all((apply(ptrans, 1, sum) >= 1 - .Machine$double.eps) | (apply(ptrans, 1, sum) <= 1 + .Machine$double.eps))) {
    stop("'ptrans' is not a stochastic matrix (column sums accross rows must be equal to one for each row)")
  }
  
  #############################
  # Checking parameter type.sojourn
  #############################
  
  type.sojourn <- match.arg(type.sojourn)
  
  #############################
  # Checking parameters distr and param
  #############################
  
  distrib.vec <- c("unif", "geom", "pois", "dweibull", "nbinom", NA)
  if (!all(distr %in% distrib.vec)) {
    stop("The specified distributions must be either ", paste(distrib.vec, collapse = ", "), 
         ".\n Incorrect distribution(s) found in 'distr': ", paste(as.character(distr)[!(distr %in% distrib.vec)], collapse = ", "))
  }
  
  if (!all(param >= 0 | is.na(param))) {
    stop("Every element of 'param' must be positive values (or NA for missing/unused parameters)")
  }
  
  if (type.sojourn == "fij") {
    
    if (!(is.matrix(distr) & dim(distr)[1] == s & dim(distr)[2] == s)) {
      stop("'distr' must be a matrix of dimension (s, s) since 'type.sojourn == \"fij\"'")
    }
    
    if (anyNA(distr[ptrans != 0])) {
      
      indexdiag <- seq(1, s * s, by = s + 1)
      distr.temp <- as.vector(distr)[-indexdiag]
      
      statesi <- row(distr)[-indexdiag][is.na(distr.temp)]
      statesj <- col(distr)[-indexdiag][is.na(distr.temp)]
      
      stop("Some conditional sojourn time distributions are not specified while 
           transitions are allowed via the matrix transition 'ptrans' (see transitions ", 
           paste0(sapply(1:length(statesi), function(x) 
             paste0("(i=", statesi[x], " to j=", statesj[x], ")")), collapse = ", "), ")")
    }
    
    if (!(all(is.na(diag(distr))))) {
      stop("All the diagonal elements of 'distr' must be equal to NA since transitions to the same state are not allowed")
    }
    
    if (!(is.array(param) & !is.matrix(param) & dim(param)[1] == s & dim(param)[2] == s)) {
      stop("'param' must be an array of dimension (s, s, 2) since 'type.sojourn == \"fij\"'")
    }
    
    if (!(all(is.na(diag(param[, , 1]))) & all(is.na(diag(param[, , 2]))))) {
      stop("All the diagonal elements of 'param' must be equal to NA since transitions to the same state are not allowed")
    }
    
    allChecking <- c()
    for (i in 1:s) {
      for (j in 1:s) {
        if (i != j & !is.na(distr[i, j])) {
          checking <- checkParameter(distr[i, j], param[i, j, ])
          if (length(checking)) {
            allChecking <- c(allChecking, paste0("-Transition (i = \"", states[i], "\" to j = \"", states[j], "\"): ", checking))
          }
        }
      }
    }
    
    if (length(allChecking)) {
      stop("Bad parameter specifications:\n\n", paste0(allChecking, collapse = "\n\n"))
    }
    
  }
  
  if (type.sojourn == "fi" | type.sojourn == "fj") {
    
    if (!(is.vector(distr) & length(distr) == s)) {
      stop("'distr' must be a vector of length s since 'type.sojourn == \"fi\"' or 'type.sojourn == \"fj\"'")
    }
    
    if (anyNA(distr)) {
      stop("'distr' cannot contain non specified distributions")
    }
    
    if (!(is.matrix(param) & dim(param)[1] == s & dim(param)[2] == 2)) {
      stop("'param' must be a matrix of dimension (s, 2) since 'type.sojourn == \"fi\"' or 'type.sojourn == \"fj\"'")
    }
    
    allChecking <- c()
    for (i in 1:s) {
      checking <- checkParameter(distr[i], param[i, ])
      if (length(checking)) {
        allChecking <- c(allChecking, paste0("-State ", ifelse(type.sojourn == "fi", "i", "j"), " = \"", states[i], "\": ", checking))
      }
    }
    
    if (length(allChecking) != 0) {
      stop("Bad parameter specifications:\n\n", paste0(allChecking, collapse = "\n\n"))
    }
    
  }
  
  if (type.sojourn == "f") {
    
    if (!(is.vector(distr) & length(distr) == 1)) {
      stop("'distr' must be one distribution since 'type.sojourn == \"f\"'")
    }
    
    if (is.na(distr)) {
      stop("'distr' must be either ", paste(distrib.vec[-length(distrib.vec)], collapse = ", "))
    }
    
    if (!(is.vector(param) & length(param) == 2)) {
      stop("'param' must be a vector of length 2 since 'type.sojourn == \"f\"'")
    }
    
    checking <- checkParameter(distr, param)
    if (length(checking) != 0) {
      stop("Bad parameter specifications :\n", checking)
    }
    
  }
  
  #############################
  # Checking parameter cens.beg and cens.end
  #############################
  
  if (!(is.logical(cens.beg) & is.logical(cens.end))) {
    stop("'cens.beg' and 'cens.end' must be TRUE or FALSE")
  }
  
  
  # Add names to the attributes init, ptrans, distr and param for readability
  colnames(ptrans) <- words(length = 1, alphabet = states)
  row.names(ptrans) <- colnames(ptrans)
  names(init) <- colnames(ptrans)
  
  if (is.matrix(distr)) {
    colnames(distr) <- colnames(ptrans)
    row.names(distr) <- colnames(distr)
    dimnames(param) <- rep(list(colnames(distr)), 2)
  } else if (length(distr) != 1) {
    names(distr) <- colnames(ptrans)
    row.names(param) <- colnames(ptrans)
  }
  
  
  ans <-
    list(
      states = states,
      s = s,
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


#' Function to check if an object is of class `smmparametric`
#' 
#' @description `is.smmparametric` returns `TRUE` if `x` is an object of 
#'   class `smmparametric`.
#' 
#' @param x An arbitrary R object.
#' @return `is.smmparametric` returns `TRUE` or `FALSE` depending on whether 
#'   `x` is an object of class `smmparametric` or not.
#' 
#' @export
#' 
is.smmparametric <- function(x) {
  inherits(x, "smmparametric")
}


# Method used to compute the semi-Markov kernel q 
# (see method getKernel.smmparametric)
.get.fijk.smmparametric <- function(x, k) {
  
  s <- x$s
  
  if (x$type.sojourn == "fij") {
    param1 <- x$param[, , 1]
    param2 <- x$param[, , 2]
    f <- matrix(0, nrow = s * s, ncol = k)
  } else if (x$type.sojourn == "fj") {
    param1 <- x$param[, 1]
    param2 <- x$param[, 2]
    f <- matrix(0, nrow = s, ncol = k)
  } else if (x$type.sojourn == "fi") {
    param1 <- x$param[, 1]
    param2 <- x$param[, 2]
    f <- matrix(0, nrow = s, ncol = k)
  } else {
    param1 <- x$param[1]
    param2 <- x$param[2]
    f <- matrix(0, nrow = 1, ncol = k)
  }
  
  if ("dweibull" %in% x$distr) {
    indices <- which(x$distr == "dweibull")
    for (j in indices) {
      f[j, ] <- ddweibull(x = 1:k, q = param1[j], beta = param2[j], zero = FALSE)
    }
  }
  if ("geom" %in% x$distr) {
    indices <- which(x$distr == "geom")
    for (j in indices) {
      f[j, ] <- dgeom(x = 0:(k - 1), prob = param1[j])
    }
  }
  if ("nbinom" %in% x$distr) {
    indices <- which(x$distr == "nbinom")
    for (j in indices) {
      f[j, ] <- dnbinom(x = 0:(k - 1), size = param1[j], prob = param2[j])
    }
  }
  if ("pois" %in% x$distr) {
    indices <- which(x$distr == "pois")
    for (j in indices) {
      f[j, ] <- dpois(x = 0:(k - 1), lambda = param1[j])
    }
  }
  if ("unif" %in% x$distr) {
    indices <- which(x$distr == "unif")
    for (j in indices) {
      f[j, ] <- sapply(1:k, function(k) ifelse(k <= x$param[j], 1 / x$param[j], 0))
    }
  }
  
  if (x$type.sojourn == "fij") {
    fijk <- array(f, c(s, s, k))
  } else if (x$type.sojourn == "fi") {
    f <- rep(as.vector(f), each = s)
    fmat <- matrix(f, nrow = k, ncol = s * s, byrow = TRUE)
    fk <- array(as.vector(t(fmat)), c(s, s, k))
    fijk <- apply(X = fk, MARGIN =  c(1, 3), FUN =  t)
  } else if (x$type.sojourn == "fj") {
    f <- rep(as.vector(f), each = s)
    fmat <- matrix(f, nrow = k, ncol = s * s, byrow = TRUE)
    fijk <- array(as.vector(t(fmat)), c(s, s, k))
  } else {
    f <- rep(f, each = s * s)
    fmat <- matrix(f, nrow = k, ncol = s * s, byrow = TRUE)
    fijk <- array(as.vector(t(fmat)), c(s, s, k))
  }
  
  return(fijk)
  
}


# Method to get the sojourn time distribution f
.get.f.smmparametric <- function(x, k) {
  
  s <- x$s
  
  if (x$type.sojourn == "fij") {
    param1 <- x$param[, , 1]
    param2 <- x$param[, , 2]
    f <- matrix(0, nrow = s * s, ncol = k)
  } else if (x$type.sojourn == "fj") {
    param1 <- x$param[, 1]
    param2 <- x$param[, 2]
    f <- matrix(0, nrow = s, ncol = k)
  } else if (x$type.sojourn == "fi") {
    param1 <- x$param[, 1]
    param2 <- x$param[, 2]
    f <- matrix(0, nrow = s, ncol = k)
  } else {
    param1 <- x$param[1]
    param2 <- x$param[2]
    f <- matrix(0, nrow = 1, ncol = k)
  }
  
  if ("dweibull" %in% x$distr) {
    indices <- which(x$distr == "dweibull")
    for (j in indices) {
      f[j, ] <- ddweibull(1:k, q = param1[j], beta = param2[j], zero = FALSE)
    }
  }
  if ("geom" %in% x$distr) {
    indices <- which(x$distr == "geom")
    for (j in indices) {
      f[j, ] <- dgeom(0:(k - 1), prob = param1[j])
    }
  }
  if ("nbinom" %in% x$distr) {
    indices <- which(x$distr == "nbinom")
    for (j in indices) {
      f[j, ] <- dnbinom(0:(k - 1), size = param1[j], prob = param2[j])
    }
  }
  if ("pois" %in% x$distr) {
    indices <- which(x$distr == "pois")
    for (j in indices) {
      f[j, ] <- dpois(0:(k - 1), lambda = param1[j])
    }
  }
  if ("unif" %in% x$distr) {
    indices <- which(x$distr == "unif")
    for (j in indices) {
      f[j, ] <- sapply(1:k, function(k) ifelse(k <= x$param[j], 1 / x$param[j], 0))
    }
  }
  
  if (x$type.sojourn == "fij") {
    f <- array(f, c(s, s, k))
  }
  
  return(f)
  
}


# Method to get the survival/reliability function Fbar
# (useful to compute the contribution to the likelihood when censoring)
.get.Fbar.smmparametric <- function(x, k) {
  
  f <- .get.f.smmparametric(x, k)
  
  if (x$type.sojourn == "fij") {
    Fbar <- 1 - apply(X = f, MARGIN = c(1, 2), cumsum)
    Fbar <- aperm(a = Fbar, perm = c(2, 3, 1))
  } else {
    Fbar <- 1 - t(apply(f, 1, cumsum))
  }
  
  return(Fbar)
  
}


# Method to get the number of parameters
# (useful for the computation of criteria such as AIC and BIC)
.getKpar.smmparametric <- function(x) {
  
  distr <- x$distr
  
  nbDweibull <- length(which(distr == "dweibull"))
  nbGeom <- length(which(distr == "geom"))
  nbNbinom <- length(which(distr == "nbinom"))
  nbPois <- length(which(distr == "pois"))
  nbUnif <- length(which(distr == "unif"))
  
  kpar <- 2 * nbDweibull + nbGeom + 2 * nbNbinom + nbPois + nbUnif
  
  return(kpar)
}


#' Log-likelihood Function
#' 
#' @description Computation of the log-likelihood for a semi-Markov model
#' 
#' @param x An object of class [smmparametric].
#' @param processes An object of class `processesSemiMarkov`.
#' 
#' @noRd
#' 
.loglik.smmparametric <- function(x, processes) {
  
  kmax <- processes$kmax
  type.sojourn <- x$type.sojourn
  cens.beg <- x$cens.beg
  cens.end <- x$cens.end
  
  #############################
  # Let's compute the log-likelihood
  #############################
  
  init <- x$init # Initial distributiob
  Nstarti <- processes$counting$Nstarti
  maskNstarti <- Nstarti != 0 & init != 0
  
  pij <- x$ptrans # Transition matrix
  Nij <- processes$counting$Nij
  maskNij <- Nij != 0 & pij != 0
  
  f <- .get.f.smmparametric(x = x, k = kmax) # Compute the sojourn time distribution
  
  
  if (type.sojourn == "fij") {
    
    Nijk <- processes$counting$Nijk
    maskNijk <- Nijk != 0 & f != 0
    
    # Uncensored log-likelihood
    loglik <- sum(Nstarti[maskNstarti] * log(init[maskNstarti])) +
      sum(Nij[maskNij] * log(pij[maskNij])) +
      sum(Nijk[maskNijk] * log(f[maskNijk]))
    
    if (cens.beg | cens.end) {# Censoring
      
      # Contribution of the first right censored time to the log-likelihood
      Fbar <- .get.Fbar.smmparametric(x = x, k = kmax)
      
      Nbijk <- processes$counting$Nbijk
      maskNbijk <- Nbijk != 0 & Fbar != 0
      
      # Contribution of the last right censored time to the log-likelihood
      Fbarj <- t(apply(X = apply(X = getKernel.smmparametric(x = x, k = kmax)[, , -1], MARGIN = c(2, 3), sum), MARGIN = 1, cumsum))
      
      Neik <- processes$counting$Neik
      maskNeik <- Neik != 0 & Fbarj != 0
      
      loglik <- loglik + (1 * cens.beg) * sum(Nbijk[maskNbijk] * log(Fbar[maskNbijk])) +
        (1 * cens.end) * sum(Neik[maskNeik] * log(Fbarj[maskNeik]))
      
    }
    
  } else if (type.sojourn == "fi") {
    
    Nik <- processes$counting$Nik
    maskNik <- Nik != 0 & f != 0
    
    # Uncensored log-likelihood
    loglik <- sum(Nstarti[maskNstarti] * log(init[maskNstarti])) +
      sum(Nij[maskNij] * log(pij[maskNij])) +
      sum(Nik[maskNik] * log(f[maskNik]))
    
    if (cens.beg | cens.end) {# Censoring
      
      # Contribution of the first and last right censored time to the log-likelihood
      Fbar <- .get.Fbar.smmparametric(x = x, k = kmax)
      
      Nbik <- processes$counting$Nbik
      maskNbik <- Nbik != 0 & Fbar != 0
      
      Neik <- processes$counting$Neik
      maskNeik <- Neik != 0 & Fbar != 0
      
      loglik <- loglik + (1 * cens.beg) * sum(Nbik[maskNbik] * log(Fbar[maskNbik])) +
        (1 * cens.end) * sum(Neik[maskNeik] * log(Fbar[maskNeik]))
      
    }
    
  } else if (type.sojourn == "fj") {
    
    Njk <- processes$counting$Njk
    maskNjk <- Njk != 0 & f != 0
    
    # Uncensored log-likelihood
    loglik <- sum(Nstarti[maskNstarti] * log(init[maskNstarti])) +
      sum(Nij[maskNij] * log(pij[maskNij])) +
      sum(Njk[maskNjk] * log(f[maskNjk]))
    
    if (cens.beg | cens.end) {
      
      # Contribution of the first right censored time to the log-likelihood
      Fbar <- .get.Fbar.smmparametric(x = x, k = kmax)
      
      Nbjk <- processes$counting$Nbjk
      maskNbjk <- Nbjk != 0 & Fbar != 0
      
      # Contribution of the last right censored time to the log-likelihood
      Fbarj <- pij %*% Fbar
      
      Neik <- processes$counting$Neik
      maskNeik <- Neik != 0 & Fbar != 0
      
      loglik <- loglik + (1 * cens.beg) * sum(Nbjk[maskNbjk] * log(Fbar[maskNbjk])) +
        (1 * cens.end) * sum(Neik[maskNeik] * log(Fbarj[maskNeik]))
      
    }
    
  } else {
    
    Nk <- processes$counting$Nk
    maskNk <- Nk != 0 & f != 0
    
    # Uncensored log-likelihood
    loglik <- sum(Nstarti[maskNstarti] * log(init[maskNstarti])) +
      sum(Nij[maskNij] * log(pij[maskNij])) +
      sum(Nk[maskNk] * log(f[maskNk]))
    
    if (cens.beg | cens.end) {# Censoring
      
      # Contribution of the first and last right censored time to the log-likelihood
      Fbar <- .get.Fbar.smmparametric(x = x, k = kmax)
      
      Nbk <- processes$counting$Nbk
      maskNbk <- Nbk != 0 & Fbar != 0
      
      Nek <- processes$counting$Nek
      maskNek <- Nek != 0 & Fbar != 0
      
      loglik <- loglik + (1 * cens.beg) * sum(Nbk[maskNbk] * log(Fbar[maskNbk])) +
        (1 * cens.end) * sum(Nek[maskNek] * log(Fbar[maskNek]))
      
    }
    
  }
  
  return(loglik)
  
}


#' Akaike Information Criterion (AIC)
#' 
#' @description Computation of the Akaike Information Criterion.
#' 
#' @param x An object of class [smmparametric].
#' @param sequences A list of vectors representing the sequences for which the 
#'   AIC will be computed based on `x`.
#' @return Value of the AIC.
#' 
#' @noRd
#' 
#' @export
#' 
aic.smmparametric <- function(x, sequences) {
  
  loglik <- loglik(x, sequences)
  
  kpar <- .getKpar(x)
  
  aic <- -2 * loglik + 2 * kpar
  
  return(aic)
  
}


#' Bayesian Information Criterion (BIC)
#' 
#' @description Computation of the Bayesian Information Criterion.
#' 
#' @param x An object of class [smmparametric].
#' @param sequences A list of vectors representing the sequences for which the 
#'   BIC will be computed based on `x`.
#' @return Value of the BIC.
#' 
#' @noRd
#' 
#' @export
#' 
bic.smmparametric <- function(x, sequences) {
  
  loglik <- loglik(x, sequences)
  
  kpar <- .getKpar(x)
  
  n <- sum(sapply(sequences, length))
  
  bic <- -2 * loglik + log(n) * kpar
  
  return(bic)
  
}


#' Method to get the semi-Markov kernel \eqn{q}
#' 
#' @description Computes the semi-Markov kernel \eqn{q_{ij}(k)}.
#' 
#' @param x An object of class [smmparametric].
#' @param k A positive integer giving the time horizon.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param klim Optional. The time horizon used to approximate the series in the
#'   computation of the mean recurrence times vector for the asymptotic 
#'   variance.
#' @return An array giving the value of \eqn{q_{ij}(k)} at each time between 0 
#'   and `k` if `var = FALSE`. If `var = TRUE`, a list containing the 
#'   following components:
#'   \itemize{
#'    \item{x: }{an array giving the value of \eqn{q_{ij}(k)} at each time 
#'      between 0 and `k`;}
#'    \item{sigma2: }{an array giving the asymptotic variance of the estimator 
#'      \eqn{\sigma_{q}^{2}(i, j, k)}.}
#'  }
#'  
#' @noRd
#' 
#' @export
#' 
getKernel.smmparametric <- function(x, k, var = FALSE, klim = 10000) {
  
  #############################
  # Checking parameters k
  #############################
  
  if (!is.numeric(k)) {
    stop("'k' must be a positive integer")
  }
  
  if ((!((k >= 0) & ((k %% 1) == 0)))) {
    stop("'k' must be a positive integer")
  }
  
  #############################
  # Checking parameters var
  #############################
  
  if (!is.logical(var)) {
    stop("'var' must be TRUE or FALSE")
  }
  
  #############################
  # Checking parameters klim
  #############################
  
  if (!is.numeric(klim)) {
    stop("'klim' must be a positive integer")
  }
  
  if ((!((klim >= 0) & ((klim %% 1) == 0)))) {
    stop("'klim' must be a positive integer")
  }
  
  
  q <- array(data = 0, dim = c(x$s, x$s, k + 1))
  
  fijk <- .get.fijk.smmparametric(x, k)
  
  if (k > 0) {
    q[, , 2:(k + 1)] <- array(x$ptrans, c(x$s, x$s, k)) * fijk  
  }
  
  if (var) {
    
    mu <- meanRecurrenceTimes(x = x, klim = klim)
    sigma2 <- array(data = mu, dim = c(x$s, x$s, k + 1)) * q * (1 - q)
    
    return(list(x = q, sigma2 = sigma2))
    
  } else {
    
    return(q)
    
  }
  
}


#' Log-likelihood Function
#' 
#' @description Computation of the log-likelihood for a semi-Markov model.
#' 
#' @param x An object of class [smmparametric].
#' @param sequences A list of vectors representing the sequences for which the 
#'   log-likelihood will be computed based on `x`.
#' @return Value of the log-likelihood.
#' 
#' @noRd
#' 
#' @export
#' 
loglik.smmparametric <- function(x, sequences) {
  
  #############################
  # Checking parameters sequences and states
  #############################
  
  if (!(is.list(sequences) & all(sapply(sequences, class) %in% c("character", "numeric")))) {
    stop("The parameter 'sequences' should be a list of vectors")
  }
  
  if (!all(unique(unlist(sequences)) %in% x$states)) {
    stop("Some states in the list of observed sequences 'sequences' are not in the state space given by the model 'x'")
  }
  
  processes <- processesSemiMarkov(sequences = sequences, states = x$states, verbose = FALSE)
  loglik <- .loglik.smmparametric(x = x, processes = processes)
  
  return(loglik)
  
}


#' @export
plot.smmparametric <- function(x, i, j, klim = NULL, ...) {
  
  #############################
  # Checking parameters i and j
  #############################
  
  if (x$type.sojourn != "f") {
    
    if (x$type.sojourn == "fi") {
      
      if (!(i %in% x$states)) {
        stop("'i' must be a state among the state space of x")
      }
      
    } else if (x$type.sojourn == "fj") {
      
      if (!(j %in% x$states)) {
        stop("'j' must be a state among the state space of x")
      }
      
    } else {
      
      if (!(i %in% x$states)) {
        stop("'i' must be a state among the state space of x")
      }
      
      if (!(j %in% x$states)) {
        stop("'j' must be a state among the state space of x")
      }
      
      if (i == j) {
        stop(paste0("The conditional distribution for the couple (i = ", i, ", j = ", j, ") doesn't exist"))
      }
      
    }
    
  }
  
  #############################
  # Checking parameter klim
  #############################
  
  if (!is.null(klim)) {
    if (!((klim > 0) & ((klim %% 1) == 0))) {
      stop("'klim' must be a strictly positive integer")
    }
  }
  
  
  if (x$type.sojourn == "fij") {
    ind.i <- which(x$states == i)
    ind.j <- which(x$states == j)
    
    param1 <- x$param[ind.i, ind.j, 1]
    param2 <- x$param[ind.i, ind.j, 2]
    dens <- x$distr[ind.i, ind.j]
    
    ylab <- bquote(f["i=" ~ .(i) ~ ", j=" ~ .(j)](k))
    main <- paste0("Sojourn time density function for the \n current state i = \"", i, "\" and the next state j = \"", j, "\"")
    
  } else if (x$type.sojourn == "fj") {
    ind.j <- which(x$states == j)
    
    param1 <- x$param[ind.j, 1]
    param2 <- x$param[ind.j, 2]
    dens <- x$distr[ind.j]
    
    ylab <- bquote(f["j=" ~ .(j)](k))
    main <- paste0("Sojourn time density function for the next state j = \"", j, "\"")
    
  } else if (x$type.sojourn == "fi") {
    ind.i <- which(x$states == i)
    
    param1 <- x$param[ind.i, 1]
    param2 <- x$param[ind.i, 2]
    dens <- x$distr[ind.i]
    
    ylab <- bquote(f["i=" ~ .(i)](k))
    main <- paste0("Sojourn time density function for the current state i = \"", i, "\"")
    
  } else {
    param1 <- x$param[1]
    param2 <- x$param[2]
    dens <- x$distr
    
    ylab <- bquote(f(k))
    main <- paste0("Sojourn time density function")
  }
  
  if (is.na(dens)) {
    stop(paste0("The conditional distribution for the couple (i = \"", i, "\", j = \"", j, "\") doesn't exist"))
  }
  
  # Compute the quantile of order alpha if klim is NULL
  alpha <- 0.95
  if (is.null(klim)) {
    klim <- do.call(what = paste0(".q", dens), args = list(alpha, param1, param2))
  }
  
  f <- do.call(what = paste0(".d", dens), args = list(1:klim, param1, param2))
  
  plot.default(x = 1:klim, y = f, xlab = "k", ylab = ylab, ...)
  title(main = main)
}


#' @export
simulate.smmparametric <- function(object, nsim = 1, seed = NULL, ...) {
  
  ###########################################################
  ###########################################################
  # The algorithm used to simulate the sequences is the following:
  # 
  # 1. Set k = 0, T_{0} = 0 and sample J_{0} from the initial distribution \alpha;
  # 2. Sample the random variable J \sim p(J_{k} , .) and set J_{k+1} = J(\omega);
  # 3. Sample the random variable X \sim F_{J_{k} J_{k+1}}(.)
  # 4. Set T_{k+1} = T_{k} + X;
  # 5. If T_{k+1} >= M, then end;
  # 6. Else, set k = k + 1 and continue to step 2.
  # 
  ###########################################################
  ###########################################################
  
  #############################
  # Checking parameter nsim
  #############################
  
  if (!all(is.numeric(nsim), is.vector(nsim), !anyNA(nsim), nsim > 0, (nsim %% 1) == 0)) {
    stop("'nsim' must be a strictly positive integer or a vector of striclty positive integers")
  }
  
  #############################
  # Checking parameter seed
  #############################
  
  if (is.null(seed)) {
    seed <- round(as.numeric(Sys.time()))
  }
  
  if (!all(is.numeric(seed), seed >= 0, (seed %% 1) == 0)) {
    stop("'seed' must be a positive integer")
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
