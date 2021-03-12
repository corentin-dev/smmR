#' Method to get the sojourn time distribution f
#' 
#' @description Computes the conditional sojourn time distribution \eqn{f(k)}, 
#'   \eqn{f_{i}(k)}, \eqn{f_{j}(k)} or \eqn{f_{ij}(k)}.
#'   
#' @param x An object of S3 class `smmfit` or `smm`.
#' @param k A positive integer giving the time horizon.
#' @return A vector, matrix or array giving the value of \eqn{f} at each time 
#'   between 0 and `k`.
#'   
#' @noRd
#' 
.get.f <- function(x, k) {
  UseMethod(".get.f", x)
}


#' Method to get the number of parameters of the semi-Markov chain
#' 
#' @description Method to get the number of parameters of the semi-Markov 
#'   chain. This method is useful for the computation of criteria such as AIC 
#'   and BIC.
#' 
#' @param x An object for which the number of parameters can be returned.
#' @return A positive integer giving the number of parameters.
#' 
#' @noRd
#' 
.getKpar <- function(x) {
  UseMethod(".getKpar", x)
}


#' Log-likelihood Function
#' 
#' @description Computation of the log-likelihood for a semi-Markov model
#' 
#' @param x An object for which the log-likelihood can be computed.
#' @param processes An object of class `processes`.
#' 
#' @noRd
#' 
.loglik <- function(x, processes) {
  UseMethod(".loglik", x)
}


#' Akaike Information Criterion (AIC)
#' 
#' @description Generic function computing the Akaike Information Criterion of 
#'   the model `x`, with the list of sequences `sequences`.
#' 
#' @param x An object for which there exists a `loglik` attribute if 
#'   `sequences = NULL` or a `loglik` method otherwise.
#' @param sequences Optional. A list of vectors representing the sequences for 
#'   which the AIC will be computed based on `x` using the method `loglik`.
#' @return Value of the AIC.
#' 
#' @export
#' 
aic <- function(x, sequences = NULL) {
  UseMethod("aic", x)
}


#' Bayesian Information Criterion (BIC)
#' 
#' @description Generic function computing the Bayesian Information Criterion 
#'   of the model `x`, with the list of sequences `sequences`.
#' 
#' @param x An object for which there exists a `loglik` attribute if 
#'   `sequences = NULL` or a `loglik` method otherwise.
#' @param sequences Optional. A list of vectors representing the sequences for 
#'   which the AIC will be computed based on `x` using the method `loglik`.
#' @return Value of the BIC.
#' 
#' @export
#' 
bic <- function(x, sequences = NULL) {
  UseMethod("bic", x)
}


#' Method to get the semi-Markov kernel \eqn{q}
#' 
#' @description Computes the semi-Markov kernel \eqn{q_{ij}(k)}.
#' 
#' @param x An object of S3 class `smmfit` or `smm`.
#' @param k A positive integer giving the time horizon.
#' @param var Logical. If `TRUE` the asymptotic variance is computed.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes] function) for the asymptotic variance.
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
#' @export
#' 
getKernel <- function(x, k, var = FALSE, klim = 10000) {
  UseMethod("getKernel", x)
}


#' Log-likelihood Function
#' 
#' @description Generic function computing the log-likelihood of the model `x`,
#'   with the list of sequences `sequences`.
#' 
#' @param x An object for which there exists a `loglik` attribute if 
#'   `sequences = NULL`. Otherwise, the log-likelihood will be computed using 
#'   the model `x` and the sequences `sequences`.
#' @param sequences Optional. A list of vectors representing the sequences for 
#'   which the log-likelihood will be computed based on `x`.
#' @return Value of the log-likelihood.
#' 
#' @export
#' 
loglik <- function(x, sequences = NULL) {
  UseMethod("loglik", x)
}


#' Mean Sojourn Times Function
#' 
#' @description The mean sojourn time is the mean time spent in each state.
#' 
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}.
#'   
#'   We are interested in investigating the mean sojourn times of a 
#'   discrete-time semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose
#'   that the evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The state of the system is given 
#'   at each instant \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}.
#'   
#'   Let \eqn{T = (T_{n})_{n \in N}} denote the successive time points when 
#'   state changes in \eqn{(Z_{n})_{n \in N}} occur and let also 
#'   \eqn{J = (J_{n})_{n \in N}} denote the successively visited states at 
#'   these time points.
#'   
#'   The mean sojourn times vector is defined as follows:
#'   
#'   \deqn{m_{i} = E[T_{1} | Z_{0} = j] = \sum_{k \geq 0} (1 - P(T_{n + 1} - T_{n} \leq k | J_{n} = j)) = \sum_{k \geq 0} (1 - H_{j}(k)),\ i \in E}
#'   
#' @param x An object of S3 class `smmfit` or `smm`.
#' @param states Vector giving the states for which the mean sojourn time 
#'   should be computed. `states` is a subset of \eqn{E}.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes] function).
#' @return A vector of length \eqn{\textrm{card}(E)} giving the values of the 
#'   mean sojourn times for each state \eqn{i \in E}.
#' 
#' @export
#' 
meanSojournTimes <- function(x, states = x$states, klim = 10000) {
  UseMethod("meanSojournTimes", x)
}


#' Method to get the mean recurrence times \eqn{\mu}
#' 
#' @description Method to get the mean recurrence times \eqn{\mu}.
#' 
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}.
#'   
#'   We are interested in investigating the mean recurrence times of a 
#'   discrete-time semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose
#'   that the evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The state of the system is given 
#'   at each instant \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}.
#'   
#'   Let \eqn{T = (T_{n})_{n \in N}} denote the successive time points when 
#'   state changes in \eqn{(Z_{n})_{n \in N}} occur and let also 
#'   \eqn{J = (J_{n})_{n \in N}} denote the successively visited states at 
#'   these time points.
#'   
#'   The mean recurrence of an arbitrary state \eqn{j \in E} is given by:
#'   
#'   \deqn{\mu_{jj} = \frac{\sum_{i \in E} \nu(i) m_{i}}{\nu(j)}}
#'   
#'   where \eqn{(\nu(1),\dots,\nu(s))} is the stationary distribution of the 
#'   embedded Markov chain \eqn{(J_{n})_{n \in N}} and \eqn{m_{i}} is the mean 
#'   sojourn time in state \eqn{i \in E} (see [meanSojournTimes] function for 
#'   the computation).
#'   
#' @param x An object of S3 class `smmfit` or `smm`.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes] function).
#' @return A vector giving the mean recurrence time 
#'   \eqn{(\mu_{i})_{i \in [1,\dots,s]}}.
#' 
#' @export
#' 
meanRecurrenceTimes <- function(x, klim = 10000) {
  UseMethod("meanRecurrenceTimes", x)
}

#' Reliability Function
#' 
#' @description Consider a system \eqn{S_{ystem}} starting to function at time 
#'   \eqn{k = 0}. The reliability or the survival function of \eqn{S_{ystem}} 
#'   at time \eqn{k \in N} is the probability that the system has functioned 
#'   without failure in the period \eqn{[0, k]}.
#'   
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}. 
#'   Denote by \eqn{U = \{1,\dots,s_1\}} the subset of operational states of 
#'   the system (the up states) and by \eqn{D = \{s_1 + 1,\dots, s\}} the 
#'   subset of failure states (the down states), with \eqn{0 < s_1 < s} 
#'   (obviously, \eqn{E = U \cup D} and \eqn{U \cap D = \emptyset}, 
#'   \eqn{U \neq \emptyset,\ D \neq \emptyset}). One can think of the states 
#'   of \eqn{U} as different operating modes or performance levels of the 
#'   system, whereas the states of \eqn{D} can be seen as failures of the 
#'   systems with different modes.
#'   
#'   We are interested in investigating the reliability of a discrete-time 
#'   semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose that the 
#'   evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The system starts to work at 
#'   instant \eqn{0} and the state of the system is given at each instant 
#'   \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}, for a certain 
#'   \eqn{i \in U}, means that the system \eqn{S_{ystem}} is in operating mode 
#'   \eqn{i} at time \eqn{k}, whereas \eqn{\{Z_k = j\}}, for a certain 
#'   \eqn{j \in D}, means that the system is not operational at time \eqn{k} 
#'   due to the mode of failure \eqn{j} or that the system is under the 
#'   repairing mode \eqn{j}.
#'   
#'   Let \eqn{T_D} denote the first passage time in subset \eqn{D}, called 
#'   the lifetime of the system, i.e.,
#'   
#'  \deqn{T_D := \textrm{inf}\{ n \in N;\ Z_n \in D\}\ \textrm{and}\ \textrm{inf}\ \emptyset := \infty.}
#'  
#'  The reliability or the survival function at time \eqn{k \in N} of a 
#'  discrete-time semi-Markov system is:
#'  
#'  \deqn{R(k) := P(T_D > k) = P(Zn \in U,n = 0,\dots,k)}
#'  
#'  which can be rewritten as follows:
#'  
#'  \deqn{R(k) = \sum_{i \in U} P(Z_0 = i) P(T_D > k | Z_0 = i) = \sum_{i \in U} \alpha_i P(T_D > k | Z_0 = i)}
#'  
#' @param x An object of S3 class `smmfit` or `smm`.
#' @param k A positive integer giving the period \eqn{[0, k]} on which the 
#'   reliability should be computed.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param level Confidence level of the asymptotic confidence interval. Helpful
#'   for an object `x` of class `smmfit`.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes][meanSojournTimes] function) for the asymptotic 
#'   variance.
#' @return A matrix with \eqn{k + 1} rows, and with columns giving values of 
#'   the reliability, variances, lower and upper asymptotic confidence limits 
#'   (if `x` is an object of class `smmfit`).
#'  
#' @references
#' V. S. Barbu, N. Limnios. (2008). Semi-Markov Chains and Hidden Semi-Markov 
#' Models Toward Applications - Their Use in Reliability and DNA Analysis. 
#' New York: Lecture Notes in Statistics, vol. 191, Springer.
#' 
#' @export
#' 
reliability <- function(x, k, upstates = x$states, level = 0.95, klim = 10000) {
  UseMethod("reliability", x)
}


#' Maintainability Function
#' 
#' @description For a reparable system \eqn{S_{ystem}} for which the failure 
#'   occurs at time \eqn{k = 0}, its maintainability at time \eqn{k \in N} is 
#'   the probability that the system is repaired up to time \eqn{k}, given that
#'   it has failed at time \eqn{k = 0}.
#'   
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}. 
#'   Denote by \eqn{U = \{1,\dots,s_1\}} the subset of operational states of 
#'   the system (the up states) and by \eqn{D = \{s_1 + 1,\dots, s\}} the 
#'   subset of failure states (the down states), with \eqn{0 < s_1 < s} 
#'   (obviously, \eqn{E = U \cup D} and \eqn{U \cap D = \emptyset}, 
#'   \eqn{U \neq \emptyset,\ D \neq \emptyset}). One can think of the states 
#'   of \eqn{U} as different operating modes or performance levels of the 
#'   system, whereas the states of \eqn{D} can be seen as failures of the 
#'   systems with different modes.
#'   
#'   We are interested in investigating the maintainability of a discrete-time 
#'   semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose that the 
#'   evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The system starts to fail at 
#'   instant \eqn{0} and the state of the system is given at each instant 
#'   \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}, for a certain 
#'   \eqn{i \in U}, means that the system \eqn{S_{ystem}} is in operating mode 
#'   \eqn{i} at time \eqn{k}, whereas \eqn{\{Z_k = j\}}, for a certain 
#'   \eqn{j \in D}, means that the system is not operational at time \eqn{k} 
#'   due to the mode of failure \eqn{j} or that the system is under the 
#'   repairing mode \eqn{j}.
#'   
#'   Thus, we take \eqn{(\alpha_{i} := P(Z_{0} = i))_{i \in U} = 0} and we 
#'   denote by \eqn{T_U} the first hitting time of subset \eqn{U}, called the 
#'   duration of repair or repair time, that is,
#'   
#'   \deqn{T_U := \textrm{inf}\{ n \in N;\ Z_n \in U\}\ \textrm{and}\ \textrm{inf}\ \emptyset := \infty.}
#'   
#'   The maintainability at time \eqn{k \in N} of a discrete-time semi-Markov 
#'   system is 
#'   
#'   \deqn{M(k) = P(T_U \leq k) = 1 - P(T_{U} > k) = 1 - P(Z_{n} \in D,\ n = 0,\dots,k).}
#'   
#' @param x An object of S3 class `smmfit` or `smm`.
#' @param k A positive integer giving the period \eqn{[0, k]} on which the 
#'   maintainability should be computed.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param level Confidence level of the asymptotic confidence interval. Helpful
#'   for an object `x` of class `smmfit`.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes] function) for the asymptotic variance.
#' @return A matrix with \eqn{k + 1} rows, and with columns giving values of 
#'   the maintainability, variances, lower and upper asymptotic confidence limits 
#'   (if `x` is an object of class `smmfit`).
#'  
#' @references
#' V. S. Barbu, N. Limnios. (2008). Semi-Markov Chains and Hidden Semi-Markov 
#' Models Toward Applications - Their Use in Reliability and DNA Analysis. 
#' New York: Lecture Notes in Statistics, vol. 191, Springer.
#' 
#' @export
#' 
maintainability <- function(x, k, upstates = x$states, level = 0.95, klim = 10000) {
  UseMethod("maintainability", x)
}


#' Availability Function
#' 
#' @description The pointwise (or instantaneous) availability of a system 
#'   \eqn{S_{ystem}} at time \eqn{k \in N} is the probability that the system 
#'   is operational at time \eqn{k} (independently of the fact that the system 
#'   has failed or not in \eqn{[0, k)}).
#' 
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}. 
#'   Denote by \eqn{U = \{1,\dots,s_1\}} the subset of operational states of 
#'   the system (the up states) and by \eqn{D = \{s_1 + 1,\dots,s\}} the 
#'   subset of failure states (the down states), with \eqn{0 < s_1 < s} 
#'   (obviously, \eqn{E = U \cup D} and \eqn{U \cap D = \emptyset}, 
#'   \eqn{U \neq \emptyset,\ D \neq \emptyset}). One can think of the states 
#'   of \eqn{U} as different operating modes or performance levels of the 
#'   system, whereas the states of \eqn{D} can be seen as failures of the 
#'   systems with different modes.
#'   
#'   We are interested in investigating the availability of a discrete-time 
#'   semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose that the 
#'   evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The state of the system is given 
#'   at each instant \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}, 
#'   for a certain \eqn{i \in U}, means that the system \eqn{S_{ystem}} is in 
#'   operating mode \eqn{i} at time \eqn{k}, whereas \eqn{\{Z_k = j\}}, for a 
#'   certain \eqn{j \in D}, means that the system is not operational at time 
#'   \eqn{k} due to the mode of failure \eqn{j} or that the system is under the
#'   repairing mode \eqn{j}.
#'   
#'   The pointwise (or instantaneous) availability of a system \eqn{S_{ystem}} 
#'   at time \eqn{k \in N} is the probability that the system is operational 
#'   at time \eqn{k} (independently of the fact that the system has failed or 
#'   not in \eqn{[0, k)}).
#'   
#'   Thus, the pointwise availability of a semi-Markov system at time 
#'   \eqn{k \in N} is
#'   
#'   \deqn{A(k) = P(Z_k \in U) = \sum_{i \in E} \alpha_i A_i(k),}
#'   
#'   where we have denoted by \eqn{A_i(k)} the conditional availability of the 
#'   system at time \eqn{k \in N}, given that it starts in state \eqn{i \in E},
#'   
#'   \deqn{A_i(k) = P(Z_k \in U | Z_0 = i).}
#'   
#' @param x An object of S3 class `smmfit` or `smm`.
#' @param k A positive integer giving the time at which the availability 
#'   should be computed.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param level Confidence level of the asymptotic confidence interval. Helpful
#'   for an object `x` of class `smmfit`.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes] function) for the asymptotic variance.
#' @return A matrix with \eqn{k + 1} rows, and with columns giving values of 
#'   the availability, variances, lower and upper asymptotic confidence limits 
#'   (if `x` is an object of class `smmfit`).
#'  
#' @references
#' V. S. Barbu, N. Limnios. (2008). Semi-Markov Chains and Hidden Semi-Markov 
#' Models Toward Applications - Their Use in Reliability and DNA Analysis. 
#' New York: Lecture Notes in Statistics, vol. 191, Springer.
#' 
#' @export
#' 
availability <- function(x, k, upstates = x$states, level = 0.95, klim = 10000) {
  UseMethod("availability", x)
}


#' BMP-Failure Rate Function
#' 
#' @description Consider a system \eqn{S_{ystem}} starting to work at time 
#'   \eqn{k = 0}. The BMP-failure rate at time \eqn{k \in N} is the conditional 
#'   probability that the failure of the system occurs at time \eqn{k}, given 
#'   that the system has worked until time \eqn{k - 1}.
#' 
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}. 
#'   Denote by \eqn{U = \{1,\dots,s_1\}} the subset of operational states of 
#'   the system (the up states) and by \eqn{D = \{s_1 + 1,\dots,s\}} the 
#'   subset of failure states (the down states), with \eqn{0 < s_1 < s} 
#'   (obviously, \eqn{E = U \cup D} and \eqn{U \cap D = \emptyset}, 
#'   \eqn{U \neq \emptyset,\ D \neq \emptyset}). One can think of the states 
#'   of \eqn{U} as different operating modes or performance levels of the 
#'   system, whereas the states of \eqn{D} can be seen as failures of the 
#'   systems with different modes.
#'   
#'   We are interested in investigating the failure rate of a discrete-time 
#'   semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose that the 
#'   evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The system starts to work at 
#'   instant \eqn{0} and the state of the system is given at each instant 
#'   \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}, for a certain 
#'   \eqn{i \in U}, means that the system \eqn{S_{ystem}} is in operating mode 
#'   \eqn{i} at time \eqn{k}, whereas \eqn{\{Z_k = j\}}, for a certain 
#'   \eqn{j \in D}, means that the system is not operational at time \eqn{k} 
#'   due to the mode of failure \eqn{j} or that the system is under the 
#'   repairing mode \eqn{j}.
#'   
#'   The BMP-failure rate at time \eqn{k \in N} is the conditional probability 
#'   that the failure of the system occurs at time \eqn{k}, given that the 
#'   system has worked until time \eqn{k - 1}.
#'   
#'   Let \eqn{T_D} denote the first passage time in subset \eqn{D}, called 
#'   the lifetime of the system, i.e.,
#'   
#'   \deqn{T_D := \textrm{inf}\{ n \in N;\ Z_n \in D\}\ \textrm{and}\ \textrm{inf}\ \emptyset := \infty.}
#'   
#'   For a discrete-time semi-Markov system, the failure rate at time 
#'   \eqn{k \geq 1} has the expression:
#'   
#'   \deqn{\lambda(k) := P(T_{D} = k | T_{D} \geq k)}
#'   
#'   We can rewrite it as follows :
#'   
#'   \deqn{\lambda(k) = 1 - \frac{R(k)}{R(k - 1)},\ \textrm{if } R(k - 1) \neq 0;\ \lambda(k) = 0, \textrm{otherwise}}
#'   
#'   The failure rate at time \eqn{k = 0} is defined by \eqn{\lambda(0) := 1 - R(0)},
#'   with \eqn{R} being the reliability function (see [reliability] function).
#'   
#'   The calculation of the reliability \eqn{R} involves the computation of 
#'   many convolutions. It implies that the computation error, may be higher 
#'   (in value) than the "true" reliability itself for reliability close to 0.
#'   In order to avoid inconsistent values of the BMP-failure rate, we use the 
#'   following formula:
#'   
#'   \deqn{\lambda(k) = 1 - \frac{R(k)}{R(k - 1)},\ \textrm{if } R(k - 1) \geq \epsilon;\ \lambda(k) = 0, \textrm{otherwise}}
#'   
#'   with \eqn{\epsilon}, the threshold, the parameter `epsilon` in the 
#'   function `failureRateBMP`.
#'   
#' @param x An object of S3 class `smmfit` or `smm`.
#' @param k A positive integer giving the period \eqn{[0, k]} on which the 
#'   BMP-failure rate should be computed.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param level Confidence level of the asymptotic confidence interval. Helpful
#'   for an object `x` of class `smmfit`.
#' @param epsilon Value of the reliability above which the latter is supposed 
#'   to be 0 because of computation errors (see Details).
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes] function) for the asymptotic variance.
#' @return A matrix with \eqn{k + 1} rows, and with columns giving values of 
#'   the BMP-failure rate, variances, lower and upper asymptotic confidence 
#'   limits (if `x` is an object of class `smmfit`).
#'  
#' @noRd
#' 
.failureRateBMP <- function(x, k, upstates = x$states, level = 0.95, epsilon = 1e-3, klim = 10000) {
  UseMethod(".failureRateBMP", x)
}


#' RG-Failure Rate Function
#' 
#' @description Discrete-time adapted failure rate, proposed by D. Roy and 
#'   R. Gupta. Classification of discrete lives. Microelectronics Reliability, 
#'   32(10):1459--1473, 1992.
#'   We call it the RG-failure rate and denote it by \eqn{r(k),\ k \in N}.
#' 
#' @details Expressing \eqn{r(k)} in terms of the reliability \eqn{R} we obtain 
#'   that the RG-failure rate function for a discrete-time system is given by:
#'   
#'   \deqn{r(k) = - \ln \frac{R(k)}{R(k - 1)},\ \textrm{if } k \geq 1;\ r(k) = - \ln R(0),\ \textrm{if } k = 0}
#'   
#'   for \eqn{R(k) \neq 0}. If \eqn{R(k) = 0}, we set \eqn{r(k) := 0}.
#'   
#'   Note that the RG-failure rate is related to the BMP-failure rate 
#'   (see [failureRateBMP] function) by:
#'   
#'   \deqn{r(k) = - \ln (1 - \lambda(k)),\ k \in N}
#'   
#'   The computation of the RG-failure rate is based on the [failureRateBMP] 
#'   function (See [failureRateBMP] for details about the parameter `epsilon`).
#'   
#' @param x An object of S3 class `smmfit` or `smm`.
#' @param k A positive integer giving the period \eqn{[0, k]} on which the 
#'   RG-failure rate should be computed.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param level Confidence level of the asymptotic confidence interval. Helpful
#'   for an object `x` of class `smmfit`.
#' @param epsilon Value of the reliability above which the latter is supposed 
#'   to be 0 because of computation errors (see Details).
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes] function) for the asymptotic variance.
#' @return A matrix with \eqn{k + 1} rows, and with columns giving values of 
#'   the RG-failure rate, variances, lower and upper asymptotic confidence 
#'   limits (if `x` is an object of class `smmfit`).
#'  
#' @noRd
#' 
.failureRateRG <- function(x, k, upstates = x$states, level = 0.95, epsilon = 1e-3, klim = 10000) {
  UseMethod(".failureRateRG", x)
}


#' Failure Rate Function
#' 
#' @description Function to compute the BMP-failure rate or the RG-failure rate.
#' 
#'   Consider a system \eqn{S_{ystem}} starting to work at time 
#'   \eqn{k = 0}. The BMP-failure rate at time \eqn{k \in N} is the conditional 
#'   probability that the failure of the system occurs at time \eqn{k}, given 
#'   that the system has worked until time \eqn{k - 1}.
#'   
#'   The RG-failure rate is a discrete-time adapted failure rate, proposed by 
#'   D. Roy and R. Gupta. Classification of discrete lives. Microelectronics 
#'   Reliability, 32(10):1459--1473, 1992. We call it the RG-failure rate and 
#'   denote it by \eqn{r(k),\ k \in N}.
#' 
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}. 
#'   Denote by \eqn{U = \{1,\dots,s_1\}} the subset of operational states of 
#'   the system (the up states) and by \eqn{D = \{s_1 + 1,\dots,s\}} the 
#'   subset of failure states (the down states), with \eqn{0 < s_1 < s} 
#'   (obviously, \eqn{E = U \cup D} and \eqn{U \cap D = \emptyset}, 
#'   \eqn{U \neq \emptyset,\ D \neq \emptyset}). One can think of the states 
#'   of \eqn{U} as different operating modes or performance levels of the 
#'   system, whereas the states of \eqn{D} can be seen as failures of the 
#'   systems with different modes.
#'   
#'   We are interested in investigating the failure rate of a discrete-time 
#'   semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose that the 
#'   evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The system starts to work at 
#'   instant \eqn{0} and the state of the system is given at each instant 
#'   \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}, for a certain 
#'   \eqn{i \in U}, means that the system \eqn{S_{ystem}} is in operating mode 
#'   \eqn{i} at time \eqn{k}, whereas \eqn{\{Z_k = j\}}, for a certain 
#'   \eqn{j \in D}, means that the system is not operational at time \eqn{k} 
#'   due to the mode of failure \eqn{j} or that the system is under the 
#'   repairing mode \eqn{j}.
#'   
#'   The BMP-failure rate at time \eqn{k \in N} is the conditional probability 
#'   that the failure of the system occurs at time \eqn{k}, given that the 
#'   system has worked until time \eqn{k - 1}.
#'   
#'   Let \eqn{T_D} denote the first passage time in subset \eqn{D}, called 
#'   the lifetime of the system, i.e.,
#'   
#'   \deqn{T_D := \textrm{inf}\{ n \in N;\ Z_n \in D\}\ \textrm{and}\ \textrm{inf}\ \emptyset := \infty.}
#'   
#'   For a discrete-time semi-Markov system, the failure rate at time 
#'   \eqn{k \geq 1} has the expression:
#'   
#'   \deqn{\lambda(k) := P(T_{D} = k | T_{D} \geq k)}
#'   
#'   We can rewrite it as follows :
#'   
#'   \deqn{\lambda(k) = 1 - \frac{R(k)}{R(k - 1)},\ \textrm{if } R(k - 1) \neq 0;\ \lambda(k) = 0, \textrm{otherwise}}
#'   
#'   The failure rate at time \eqn{k = 0} is defined by \eqn{\lambda(0) := 1 - R(0)},
#'   with \eqn{R} being the reliability function (see [reliability] function).
#'   
#'   The calculation of the reliability \eqn{R} involves the computation of 
#'   many convolutions. It implies that the computation error, may be higher 
#'   (in value) than the "true" reliability itself for reliability close to 0.
#'   In order to avoid inconsistent values of the BMP-failure rate, we use the 
#'   following formula:
#'   
#'   \deqn{\lambda(k) = 1 - \frac{R(k)}{R(k - 1)},\ \textrm{if } R(k - 1) \geq \epsilon;\ \lambda(k) = 0, \textrm{otherwise}}
#'   
#'   with \eqn{\epsilon}, the threshold, the parameter `epsilon` in the 
#'   function `failureRate`.
#'   
#'   
#'   Expressing the RG-failure rate \eqn{r(k)} in terms of the reliability 
#'   \eqn{R} we obtain that the RG-failure rate function for a discrete-time 
#'   system is given by:
#'   
#'   \deqn{r(k) = - \ln \frac{R(k)}{R(k - 1)},\ \textrm{if } k \geq 1;\ r(k) = - \ln R(0),\ \textrm{if } k = 0}
#'   
#'   for \eqn{R(k) \neq 0}. If \eqn{R(k) = 0}, we set \eqn{r(k) := 0}.
#'   
#'   Note that the RG-failure rate is related to the BMP-failure rate by:
#'   
#'   \deqn{r(k) = - \ln (1 - \lambda(k)),\ k \in N}
#'   
#' @param x An object of S3 class `smmfit` or `smm`.
#' @param k A positive integer giving the period \eqn{[0, k]} on which the 
#'   BMP-failure rate should be computed.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param failure.rate Type of failure rate to compute. If `failure.rate = "BMP"` 
#'   (default value), then BMP-failure-rate is computed. If `failure.rate = "RG"`, 
#'   the RG-failure rate is computed.
#' @param level Confidence level of the asymptotic confidence interval. Helpful
#'   for an object `x` of class `smmfit`.
#' @param epsilon Value of the reliability above which the latter is supposed 
#'   to be 0 because of computation errors (see Details).
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes] function) for the asymptotic variance.
#' @return A matrix with \eqn{k + 1} rows, and with columns giving values of 
#'   the BMP-failure rate or RG-failure rate, variances, lower and upper 
#'   asymptotic confidence limits (if `x` is an object of class `smmfit`).
#'  
#' @references
#' V. S. Barbu, N. Limnios. (2008). Semi-Markov Chains and Hidden Semi-Markov 
#' Models Toward Applications - Their Use in Reliability and DNA Analysis. 
#' New York: Lecture Notes in Statistics, vol. 191, Springer.
#' 
#' R.E. Barlow, A.W. Marshall, and F. Prochan. (1963). Properties of probability 
#' distributions with monotone hazard rate. Ann. Math. Statist., 34, 375-389.
#' 
#' D. Roy and R. Gupta. (1992). Classification of discrete lives. 
#' Microelectron. Reliabil., 32 (10), 1459-1473.
#' 
#' @export
#' 
failureRate <- function(x, k, upstates = x$states, failure.rate = c("BMP", "RG"), level = 0.95, epsilon = 1e-3, klim = 10000) {
  UseMethod("failureRate", x)
}


#' Mean Time To Failure (MTTF) Function
#' 
#' @description Consider a system \eqn{S_{ystem}} starting to work at time 
#'   \eqn{k = 0}. The mean time to failure (MTTF) is defined as the mean 
#'   lifetime.
#' 
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}. 
#'   Denote by \eqn{U = \{1,\dots,s_1\}} the subset of operational states of 
#'   the system (the up states) and by \eqn{D = \{s_1 + 1,\dots,s\}} the 
#'   subset of failure states (the down states), with \eqn{0 < s_1 < s} 
#'   (obviously, \eqn{E = U \cup D} and \eqn{U \cap D = \emptyset}, 
#'   \eqn{U \neq \emptyset,\ D \neq \emptyset}). One can think of the states 
#'   of \eqn{U} as different operating modes or performance levels of the 
#'   system, whereas the states of \eqn{D} can be seen as failures of the 
#'   systems with different modes.
#'   
#'   We are interested in investigating the mean time to failure of a 
#'   discrete-time semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose
#'   that the evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The system starts to work at 
#'   instant \eqn{0} and the state of the system is given at each instant 
#'   \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}, for a certain 
#'   \eqn{i \in U}, means that the system \eqn{S_{ystem}} is in operating mode 
#'   \eqn{i} at time \eqn{k}, whereas \eqn{\{Z_k = j\}}, for a certain 
#'   \eqn{j \in D}, means that the system is not operational at time \eqn{k} 
#'   due to the mode of failure \eqn{j} or that the system is under the 
#'   repairing mode \eqn{j}.
#'   
#'   Let \eqn{T_D} denote the first passage time in subset \eqn{D}, called 
#'   the lifetime of the system, i.e.,
#'   
#'   \deqn{T_D := \textrm{inf}\{ n \in N;\ Z_n \in D\}\ \textrm{and}\ \textrm{inf}\ \emptyset := \infty.}
#'   
#'   The mean time to failure (MTTF) is defined as the mean lifetime, i.e., the
#'   expectation of the hitting time to down set \eqn{D},
#'   
#'   \deqn{MTTF = E[T_{D}]}
#'   
#' @param x An object of S3 class `smmfit` or `smm`.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param level Confidence level of the asymptotic confidence interval. Helpful
#'   for an object `x` of class `smmfit`.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes] function) for the asymptotic variance.
#' @return A matrix with \eqn{\textrm{card}(U) = s_{1}} rows, and with columns 
#'   giving values of the mean time to failure for each state \eqn{i \in U}, 
#'   variances, lower and upper asymptotic confidence limits (if `x` is an 
#'   object of class `smmfit`).
#' 
#' @references
#' V. S. Barbu, N. Limnios. (2008). Semi-Markov Chains and Hidden Semi-Markov 
#' Models Toward Applications - Their Use in Reliability and DNA Analysis. 
#' New York: Lecture Notes in Statistics, vol. 191, Springer.
#' 
#' I. Votsi & A. Brouste (2019) Confidence interval for the mean time to 
#' failure in semi-Markov models: an application to wind energy production, 
#' Journal of Applied Statistics, 46:10, 1756-1773
#' 
#' @export
#' 
mttf <- function(x, upstates = x$states, level = 0.95, klim = 10000) {
  UseMethod("mttf", x)
}


#' Mean Time To Repair (MTTR) Function
#' 
#' @description Consider a system \eqn{S_{ystem}} that has just failed at time 
#'   \eqn{k = 0}. The mean time to repair (MTTR) is defined as the mean of the 
#'   repair duration.
#'   
#' @details Consider a system (or a component) \eqn{S_{ystem}} whose possible 
#'   states during its evolution in time are \eqn{E = \{1,\dots,s\}}. 
#'   Denote by \eqn{U = \{1,\dots,s_1\}} the subset of operational states of 
#'   the system (the up states) and by \eqn{D = \{s_1 + 1,\dots,s\}} the 
#'   subset of failure states (the down states), with \eqn{0 < s_1 < s} 
#'   (obviously, \eqn{E = U \cup D} and \eqn{U \cap D = \emptyset}, 
#'   \eqn{U \neq \emptyset,\ D \neq \emptyset}). One can think of the states 
#'   of \eqn{U} as different operating modes or performance levels of the 
#'   system, whereas the states of \eqn{D} can be seen as failures of the 
#'   systems with different modes.
#'   
#'   We are interested in investigating the mean time to repair of a 
#'   discrete-time semi-Markov system \eqn{S_{ystem}}. Consequently, we suppose
#'   that the evolution in time of the system is governed by an E-state space 
#'   semi-Markov chain \eqn{(Z_k)_{k \in N}}. The system has just failed at 
#'   instant 0 and the state of the system is given at each instant 
#'   \eqn{k \in N} by \eqn{Z_k}: the event \eqn{\{Z_k = i\}}, for a certain 
#'   \eqn{i \in U}, means that the system \eqn{S_{ystem}} is in operating mode 
#'   \eqn{i} at time \eqn{k}, whereas \eqn{\{Z_k = j\}}, for a certain 
#'   \eqn{j \in D}, means that the system is not operational at time \eqn{k} 
#'   due to the mode of failure \eqn{j} or that the system is under the 
#'   repairing mode \eqn{j}.
#'   
#'   Let \eqn{T_U} denote the first passage time in subset \eqn{U}, called the 
#'   duration of repair or repair time, i.e.,
#'   
#'   \deqn{T_U := \textrm{inf}\{ n \in N;\ Z_n \in U\}\ \textrm{and}\ \textrm{inf}\ \emptyset := \infty.}
#'   
#'   The mean time to repair (MTTR) is defined as the mean of the repair 
#'   duration, i.e., the expectation of the hitting time to up set \eqn{U},
#'   
#'   \deqn{MTTR = E[T_{U}]}
#'   
#' @param x An object of S3 class `smmfit` or `smm`.
#' @param upstates Vector giving the subset of operational states \eqn{U}.
#' @param level Confidence level of the asymptotic confidence interval. Helpful
#'   for an object `x` of class `smmfit`.
#' @param klim Optional. The time horizon used to approximate the series in the 
#'   computation of the mean sojourn times vector \eqn{m} (cf. 
#'   [meanSojournTimes] function) for the asymptotic variance.
#' @return A matrix with \eqn{\textrm{card}(U) = s_{1}} rows, and with columns 
#'   giving values of the mean time to repair for each state \eqn{i \in U}, 
#'   variances, lower and upper asymptotic confidence limits (if `x` is an 
#'   object of class `smmfit`).
#'  
#' @references
#' V. S. Barbu, N. Limnios. (2008). Semi-Markov Chains and Hidden Semi-Markov 
#' Models Toward Applications - Their Use in Reliability and DNA Analysis. 
#' New York: Lecture Notes in Statistics, vol. 191, Springer.
#' 
#' I. Votsi & A. Brouste (2019) Confidence interval for the mean time to 
#' failure in semi-Markov models: an application to wind energy production, 
#' Journal of Applied Statistics, 46:10, 1756-1773
#' 
#' @export
#' 
mttr <- function(x, upstates = x$states, level = 0.95, klim = 10000) {
  UseMethod("mttr", x)
}
