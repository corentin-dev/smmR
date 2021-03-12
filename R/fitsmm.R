#' Maximum Likelihood Estimation (MLE) of a semi-Markov chain
#' 
#' @description Maximum Likelihood Estimation of a semi-Markov chain starting 
#'   from one or several sequences. This estimation can be parametric or 
#'   non-parametric, non-censored, censored at the beginning and/or at the end 
#'   of the sequence, with one or several trajectories. Several parametric 
#'   distributions are considered (Uniform, Geometric, Poisson, Discrete 
#'   Weibull and Negative Binomial).
#'   
#' @details This function estimates a semi-Markov model in parametric or 
#' non-parametric case, taking into account the type of sojourn time and the 
#' censoring described in references. The non-parametric estimation concerns 
#' sojourn time distributions defined by the user. For the parametric 
#' estimation, several discrete distributions are considered (see below).
#' 
#' The difference between the Markov model and the semi-Markov model concerns 
#' the modeling of the sojourn time. With a Markov chain, the sojourn time 
#' distribution is modeled by a Geometric distribution (in discrete time). 
#' With a semi-Markov chain, the sojourn time can be any arbitrary distribution.
#' In this package, the available distribution for a semi-Markov model are :
#'  \itemize{
#'    \item Uniform: \eqn{f(x) = \frac{1}{n}} for \eqn{1 \le x \le n}. \eqn{n} is the parameter;
#'    \item Geometric: \eqn{f(x) = \theta (1-\theta)^{x - 1}} for \eqn{x = 1, 2,\ldots,n}, \eqn{0 < \theta < 1}, \eqn{\theta} is the probability of success.
#'      \eqn{\theta} is the parameter;
#'    \item Poisson: \eqn{f(x) = \frac{\lambda^x exp(-\lambda)}{x!}} for \eqn{x = 0, 1, 2,\ldots,n}, with \eqn{\lambda > 0}.
#'      \eqn{\lambda} is the parameter;
#'    \item Discrete Weibull of type 1: \eqn{f(x)=q^{(x-1)^{\beta}}-q^{x^{\beta}}}, \eqn{x = 1, 2,\ldots,n}, 
#'      with \eqn{0 < q < 1}, the first parameter and \eqn{\beta > 0} the second parameter.
#'      \eqn{(q, \beta)} are the parameters;
#'    \item Negative binomial: \eqn{f(x)=\frac{\Gamma(x+\alpha)}{\Gamma(\alpha) x!} p^{\alpha} (1 - p)^x}, 
#'      for \eqn{x = 0, 1, 2,\ldots,n}, \eqn{\Gamma} is the Gamma function, 
#'      \eqn{\alpha} is the parameter of overdispersion and \eqn{p} is the 
#'      probability of success, \eqn{0 < p < 1}. \eqn{(\alpha, p)} are the parameters;
#'    \item Non-parametric.
#'  }
#'  
#'  
#' We define :
#'  \itemize{
#'    \item the semi-Markov kernel \eqn{q_{ij}(k) = P( J_{m+1} = j, T_{m+1} - T_{m} = k | J_{m} = i )};
#'    \item the transition matrix \eqn{(p_{trans}(i,j))_{i,j} \in states} of the 
#'      embedded Markov chain \eqn{J = (J_m)_m}, \eqn{p_{trans}(i,j) = P( J_{m+1} = j | J_m = i )};
#'    \item the initial distribution \eqn{\mu_i = P(J_1 = i) = P(Z_1 = i)}, \eqn{i \in 1, 2, \dots, s};
#'    \item the conditional sojourn time distributions \eqn{(f_{ij}(k))_{i,j} \in states,\ k \in N ,\ f_{ij}(k) = P(T_{m+1} - T_m = k | J_m = i, J_{m+1} = j )},
#'      \eqn{f} is specified by the argument `param` in the parametric case 
#'      and by `distr` in the non-parametric case.
#'  }
#'  
#' The maximum likelihood estimation of the transition matrix of the embedded 
#' Markov chain is \eqn{\widehat{p_{trans}}(i,j) = \frac{N_{ij}}{N_{i.}}}.
#' 
#' Five methods are proposed for the estimation of the initial distribution :
#' \describe{
#'   \item{Maximum Likelihood Estimator: }{The Maximum Likelihood Estimator 
#'     for the initial distribution. The formula is: 
#'     \eqn{\widehat{\mu_i} = \frac{Nstart_i}{L}}, where \eqn{Nstart_i} is 
#'     the number of occurences of the word \eqn{i} (of length \eqn{k}) at 
#'     the beginning of each sequence and \eqn{L} is the number of sequences. 
#'     This estimator is reliable when the number of sequences \eqn{L} is high.}
#'   \item{Limit (stationary) distribution: }{The limit (stationary) 
#'     distribution of the semi-Markov chain is used as a surrogate of the 
#'     initial distribution.}
#'   \item{Frequencies of each state: }{The initial distribution is replaced 
#'     by taking the frequencies of each state in the sequences.}
#'   \item{Uniform distribution: }{The initial probability of each state is 
#'     equal to \eqn{1 / s}, with \eqn{s}, the number of states.}
#'  }
#'  
#' Note that \eqn{q_{ij}(k) = p_{trans}(i,j) \ f_{ij}(k)} in the general case 
#' (depending on the present state and on the next state). For particular cases, 
#' we replace \eqn{f_{ij}(k)} by \eqn{f_{i.}(k)} (depending on the present 
#' state \eqn{i}), \eqn{f_{.j}(k)} (depending on the next state \eqn{j}) and 
#' \eqn{f_{..}(k)} (depending neither on the present state nor on the next 
#' state).
#' 
#' In this package we can choose different types of sojourn time. 
#' Four options are available for the sojourn times:
#' \itemize{
#'   \item depending on the present state and on the next state (`fij`);
#'   \item depending only on the present state (`fi`);
#'   \item depending only on the next state (`fj`);
#'   \item depending neither on the current, nor on the next state (`f`).
#' }
#' 
#' If  `type.sojourn = "fij"`, `distr` is a matrix of dimension \eqn{(s, s)} 
#' (e.g., if the 1st element of the 2nd column is `"pois"`, that is to say we 
#' go from the first state to the second state following a Poisson distribution).
#' If `type.sojourn = "fi"` or `"fj"`, `distr` must be a vector (e.g., if the 
#' first element of vector is `"geom"`, that is to say we go from (or to) the 
#' first state to (or from) any state according to a Geometric distribution).
#' If `type.sojourn = "f"`, `distr` must be one of `"unif"`, `"geom"`, `"pois"`, 
#' `"dweibull"`, `"nbinom"` (e.g., if `distr` is equal to `"nbinom"`, that is 
#' to say that the sojourn time when going from one state to another state 
#' follows a Negative Binomial distribution).
#' For the non-parametric case, `distr` is equal to `"nonparametric"` whatever 
#' type of sojourn time given.
#' 
#' If the sequence is censored at the beginning and/or at the end, `cens.beg` 
#' must be equal to `TRUE` and/or `cens.end` must be equal to `TRUE`. 
#' All the sequences must be censored in the same way.
#' 
#' @param sequences A list of vectors representing the sequences.
#' @param states Vector of state space (of length s).
#' @param type.sojourn Type of sojourn time (for more explanations, see Details).
#' @param distr By default `"nonparametric"` for the non-parametric estimation 
#'   case.
#'   
#'   If the parametric estimation case is desired, `distr` should be:
#'   \itemize{
#'     \item A matrix of distributions of dimension \eqn{(s, s)} if `type.sojourn = "fij"`;
#'     \item A vector of distributions of length \eqn{s} if `type.sojourn = "fi"` or `"fj"`;
#'     \item A distribution if `type.sojourn = "f"`.
#'   }
#'   
#'   The distributions to be used in `distr` must be one of `"unif"`, `"geom"`, 
#'   `"pois"`, `"dweibull"`, `"nbinom"`.
#' @param init.estim Optional. `init.estim` gives the method used to estimate 
#'   the initial distribution. The following methods are proposed:
#'   \itemize{
#'     \item `init.estim = "mle"`: the classical Maximum Likelihood Estimator 
#'       is used to estimate the initial distribution `init`;
#'     \item `init.estim = "limit"`: the initial distribution is replaced by 
#'       the limit (stationary) distribution of the semi-Markov chain;
#'     \item `init.estim = "freq"`: the initial distribution is replaced by 
#'       the frequencies of each state in the sequences;
#'     \item `init.estim = "unif"`: the initial probability of each state is 
#'       equal to \eqn{1 / s}, with \eqn{s} the number of states.
#'   }
#' @param cens.beg Optional. A logical value indicating whether or not 
#'   sequences are censored at the beginning.
#' @param cens.end Optional. A logical value indicating whether or not 
#'   sequences are censored at the end.
#' @return Returns an object of S3 class `smmfit` (inheriting from the S3 
#'   class `smm` and [smmnonparametric] class if `distr = "nonparametric"`
#'   or [smmparametric] otherwise).
#'   
#' @seealso [smmnonparametric], [smmparametric], [simulate.smm],
#'   [simulate.smmfit], [plot.smm], [plot.smmfit]
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
#' 
#' seq1 <- simulate(object = smm1, nsim = c(1000, 10000, 2000), seed = 100)
#' 
#' # Estimation of simulated sequences
#' est1 <- fitsmm(sequences = seq1, states = states, type.sojourn = "fij", 
#'                            distr = distr.matrix)
#' est1
#' 
fitsmm <- function(sequences, states, type.sojourn = c("fij", "fi", "fj", "f"), 
                               distr = "nonparametric", init.estim = "mle", cens.beg = FALSE, cens.end = FALSE) {
  
  
  #############################
  # Checking parameters sequences and states
  #############################
  
  if (!(is.list(sequences) & all(sapply(sequences, class) %in% c("character", "numeric")))) {
    stop("The parameter 'sequences' should be a list of vectors")
  }
  
  if (!all(unique(unlist(sequences)) %in% states)) {
    stop("Some states in the list of observed sequences 'sequences' are not in the state space 'states'")
  }
  
  #############################
  # Checking parameter type.sojourn
  #############################
  
  type.sojourn <- match.arg(type.sojourn)
  
  #############################
  # Checking parameters distr
  #############################
  
  s <- length(states) # State space size
  
  if (!all(is.vector(distr), length(distr) == 1, distr == "nonparametric")) {
    
    distrib.vec <- c("unif", "geom", "pois", "dweibull", "nbinom", NA)
    if (!all(distr %in% distrib.vec)) {
      stop("The specified distributions must be either ", paste(distrib.vec, collapse = ", "), 
           ".\n Incorrect distribution(s) found in 'distr': ", paste(as.character(distr)[!(distr %in% distrib.vec)], collapse = ", "))
    }
    
    if (type.sojourn == "fij") {
      
      if (!(is.matrix(distr) & dim(distr)[1] == s & dim(distr)[2] == s)) {
        stop("'distr' must be a matrix of dimension (s, s) since 'type.sojourn == \"fij\"'")
      }
      
      if (!(all(is.na(diag(distr))))) {
        stop("All the diagonal elements of 'distr' must be equal to NA since transitions to the same state are not allowed")
      }
      
    }
    
    if (type.sojourn == "fi" | type.sojourn == "fj") {
      
      if (!(is.vector(distr) & length(distr) == s)) {
        stop("'distr' must be a vector of length s since 'type.sojourn == \"fi\"' or 'type.sojourn == \"fj\"'")
      }
      
      if (anyNA(distr)) {
        stop("'distr' cannot contain non specified distributions")
      }
      
    }
    
    if (type.sojourn == "f") {
      
      if (!(is.vector(distr) & length(distr) == 1)) {
        stop("'distr' must be one distribution since 'type.sojourn == \"f\"'")
      }
      
      if (is.na(distr)) {
        stop("'distr' must be either ", paste(distrib.vec[-length(distrib.vec)], collapse = ", "))
      }
      
    }
    
  }
  
  #############################
  # Checking parameter init.estim
  #############################
  
  # init.estim <- match.arg(init.estim)
  
  
  processes <- processesSemiMarkov(sequences = sequences, states = states, verbose = TRUE)
  
  if (all(length(distr) == 1, distr == "nonparametric")) {
    
    if (!(cens.beg)) {
      
      if (processes$L > 1) {# If more than one sequence, use the estimation based on a couple Markov chain
        smm <- .fit.nonparam.couplemarkovchain(processes = processes, states = states, type.sojourn = type.sojourn, init.estim = init.estim, cens.end = cens.end)
      } else {
        if (!cens.end) {
          smm <- .fit.nonparam.nocensoring(processes = processes, type.sojourn = type.sojourn, init.estim = init.estim, cens.beg = cens.beg)
        } else {
          smm <- .fit.nonparam.couplemarkovchain(processes = processes, states = states, type.sojourn = type.sojourn, init.estim = init.estim, cens.end = cens.end)
        }
      }
      
    } else {
      
      if (!cens.end) {
        stop("fitsmm not implemented in the case 'distr = \"nonparametric\"', 'cens.beg = TRUE', 'cens.end = FALSE'")
      } else {
        stop("fitsmm not implemented in the case 'distr = \"nonparametric\"', 'cens.beg = TRUE', 'cens.end = TRUE'")
      }
      
    }
  } else {
    smm <- .fit.param(processes = processes, states = states, type.sojourn = type.sojourn, distr = distr, init.estim = init.estim, cens.end = cens.end, cens.beg = cens.beg)
  }
  
  if (any(smm$init == 0)) {
    message("The probabilities of the initial state(s) \"", 
            paste0(names(which(smm$init == 0)), collapse = "\", \""),
            "\" are 0.")
  }
  
  loglik <- .loglik(x = smm, processes = processes)
  estimate <- smmfit(smm = smm, M = processes$M, loglik = loglik, sequences = sequences)
  
  return(estimate)
  
}
