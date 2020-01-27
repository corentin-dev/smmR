\name{simulSM}
\alias{simulSM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of a semi-Markov chain
}
\description{
Simulation of a semi-Markov chain starting from chosen parameters. This simulation can be parametric or non-parametric, non-censored, censored at the beginning and/or at the end of the sequence, with one or several trajectories. Several parametric distributions are considered (Uniform, Geometric, Poisson, Discrete Weibull and Negative Binomial). 
}
\usage{
# parametric case  
simulSM(E, NbSeq, lengthSeq, TypeSojournTime = "fij", init, Ptrans, distr, param, 
laws = NULL, cens.beg = 0, cens.end = 0, File.out = NULL)

# non-parametric case
simulSM(E, NbSeq, lengthSeq, TypeSojournTime = "fij", init, Ptrans, distr = "NP", 
param = NULL, laws, cens.beg = 0, cens.end = 0, File.out = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{
Vector of state space
}
  \item{NbSeq}{
Number of simulated sequences  
}
  \item{lengthSeq}{
Vector containing the lengths of each simulated sequence  (the length of sequences can be different, see details)
}
  \item{TypeSojournTime}{
Character : must be one of "fij", "fi", "fj", "f" (for more explanations, see Details)
}
  \item{init}{
Vector of initial distribution of length S, with S = length(E)
}
  \item{Ptrans}{
Matrix of transition probabilities of the embedded Markov chain \eqn{J=(J_m)_{m}} (The sum of the probabilities on the same row must be equal to 1. According to semi-Markov model, the probability to stay in the same state is assumed to be equal to 0.)
}
     \item{distr}{
- "NP" for nonparametric case, laws have to be used, param is useless

- Matrix of distributions of size SxS if TypeSojournTime is equal to "fij";

- Vector of distributions of size S if TypeSojournTime is equal to "fi" or "fj";

- A distribution if TypeSojournTime is equal to "f".

The distributions to be used in distr must be one of "uniform", "geom", "pois", "dweibull", "nbinom".
}

  \item{param}{
- Useless if distr = "NP"

- Array of distribution parameters of size SxSx2 (2 corresponds to the maximal number of distribution parameters) if TypeSojournTime is equal to "fij";

- Matrix of distribution parameters of size Sx2 if TypeSojournTime is equal to "fi" or "fj";

- Vector of distribution parameters of length 2 if TypeSojournTime is equal to "f". 
}
  \item{laws}{
- Useless if distr \eqn{\neq} "NP"

- Array of size SxSxKmax if TypeSojournTime is equal to "fij";

- Matrix of size SxKmax if TypeSojournTime is equal to "fi" or "fj";

- Vector of length Kmax if the TypeSojournTime is equal to "f". 

Kmax is the maximum length of the sojourn times. 
}
  \item{cens.beg}{
Type of censoring at the beginning of sample paths;  1 (if the first sojourn time is censored) or 0 (if not). All the sequences must be censored in
the same way.}
  \item{cens.end}{
Type of censoring at the end of sample paths;  1 (if the last sojourn time is censored) or 0 (if not). All the sequences must be censored in
the same way.
}
  \item{File.out}{
Name of fasta file to save the sequences;  if File.out = NULL, no file is created.
}
}
\details{
This function simulates a semi-Markov model in the parametric and non-parametric case, 
taking into account the type of sojourn time and the censoring described in references.

The non-parametric simulation concerns sojourn time distributions defined by the user.
For the parametric simulation, several discrete distributions are considered (see below).

The difference between the Markov model and the semi-Markov model concerns the modelling of the sojourn time. 
With a Markov chain, the sojourn time distribution is modeled by a Geometric distribution. 
With a semi-Markov chain, the sojourn time can be arbitrarily  distributed.
In this package, the available distribution for a semi-Markov model are:
\itemize{
\item Uniform: \eqn{f(x) = 1/n} for \eqn{a \le x \le b}, with \eqn{n = b-a+1}
\item Geometric: \eqn{f(x) = \theta (1-\theta)^x} for \eqn{x = 0, 1, 2,\ldots,n}, \eqn{0 < \theta \le 1}, with \eqn{n > 0} and \eqn{\theta} is the probability of success
\item Poisson: \eqn{f(x) = (\lambda^x exp(-\lambda))/x!} for \eqn{x = 0, 1, 2,\ldots,n}, with \eqn{n > 0} and \eqn{\lambda > 0} 
\item Discrete Weibull of type 1: \eqn{f(x)=q^{(x-1)^{\beta}}-q^{x^{\beta}}, x=1,2,3,\ldots,n}, with \eqn{n > 0}, \eqn{q} is the first parameter and \eqn{\beta} is the second parameter
\item Negative Binomial: \eqn{f(x)=\Gamma(x+\alpha)/(\Gamma(\alpha) x!) (\alpha/(\alpha+\mu))^{\alpha} (\mu/(\alpha+\mu))^x}, for \eqn{x = 0, 1, 2,\ldots,n}, \eqn{n > 0,\ \Gamma } is the gamma function, \eqn{\alpha} is the parameter of overdispersion and \eqn{\mu} is the mean
\item Non-parametric
}
We define :
\itemize{
\item the semi-Markov kernel \eqn{q_{ij}(k) = P( J_{m+1} = j, T_{m+1} - T_{m} = k | J_{m} = i )};
\item the transition matrix \eqn{(Ptrans(i,j))_{i,j} \in E} of the embedded Markov chain \eqn{J = (J_m)_m}, \eqn{Ptrans(i,j) = P( J_{m+1} = j | J_m = i )};
\item the initial distribution, \eqn{init = P(J_1 = i) = P(Y_1 = i)};
\item the conditional sojourn time distributions \eqn{(f_{ij}(k))_{i,j} \in E,\ k \in N ,\ f_{ij}(k) = P(T_{m+1} - T_m = k | J_m = i, J_{m+1} = j )}, f is specified by the argument "param" in the parametric case and by "laws" in the non-parametric case.
}
In this package we can choose differents types of sojourn times. Four options are available for the sojourn times
\itemize{
\item depending on the present state and on the next state ("fij");
\item depending only on the present state ("fi");
\item depending only on the next state ("fj");
\item depending neither on the current, nor on the next state ("f").
}

If  TypeSojournTime is equal to "fij", distr is a matrix (SxS) (e.g., if the row 1 of the 2nd column is "pois", that is to say we go from the first state to the second state following a Poisson distribution). 
If TypeSojournTime is equal to "fi" or "fj", distr must be a vector (e.g., if the first element of vector is "geom", that is to say we go from the first state to any state according to a Geometric distribution). 
If TypeSojournTime is equal to "f", distr must be one of "uniform", "geom", "pois", "dweibull", "nbinom" (e.g., if distr is equal to "nbinom", that is to say that the sojourn times when going from any state to any state follows a Negative Binomial distribution).
For the non-parametric case, distr is equal to "NP" whatever type of sojourn time given.

If the sequence is censored at the beginning and at the end, cens.beg must be equal to 1 and cens.end must be equal to 1 too. All the sequences must be censored in the same way.

Moreover, in the non-parametric case TypeSojournTime is equal to "fij" then laws must be an array of size SxSxKmax. 
If distr is equal to "NP" and TypeSojournTime is equal to "fi" or "fj" then laws must be a matrix of size SxKmax. 
If the distr is equal to "NP" and TypeSojournTime is equal to "f" then laws is a vector of length Kmax.

For the simulation of a non-censored sequence, the length of the sequence can be greater than lengthSeq.
}
\value{
simulSM returns sequences of size lengthSeq or greater

If File.out is not NULL, a fasta file will be created containing the sequences.
}

\author{
Vlad Stefan Barbu, barbu@univ-rouen.fr \cr Caroline Berard, caroline.berard@univ-rouen.fr \cr Dominique Cellier, dominique.cellier@laposte.net \cr Mathilde Sautreuil, mathilde.sautreuil@etu.univ-rouen.fr \cr Nicolas Vergne, nicolas.vergne@univ-rouen.fr 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
estimSM, simulMk
}
\examples{
alphabet = c("a","c","g","t")
S = length(alphabet)
## creation of the initial distribution
vect.init = c(1/4,1/4,1/4,1/4)
# creation of the transition matrix
Pij = matrix(c(0,0.2,0.3,0.5,0.4,0,0.2,0.4,0.1,0.2,0,0.7,0.8,0.1,0.1,0), nrow = S, 
ncol = S, byrow = TRUE)

################################
## Parametric simulation of a trajectory (of length equal to 50), 
## where the sojourn times depend neither on the present state nor on the next state ("f").
################################
## Simulation of a sequence of length 50
seq50 = simulSM(E = alphabet, NbSeq = 1, lengthSeq = 50, TypeSojournTime = "f", 
                init = vect.init, Ptrans = Pij, distr = "pois", param = 2)
          
\donttest{      
################################
## Parametric simulation of several trajectories (3 trajectories of length 1000, 
## 10 000 and  2000 respectively),
## where the sojourn times depend on the present state and on the next state ("fij")
## the sojourn time is modelled by different distributions
################################
lengthSeq3 = c(1000, 10000, 2000)
## creation of the distribution matrix
distr.matrix = matrix(c("dweibull", "pois", "geom", "nbinom", "geom", "nbinom", 
"pois", "dweibull", "pois", "pois", "dweibull", "geom", "pois","geom", "geom", 
"nbinom"), nrow = S, ncol = S, byrow = TRUE)
## creation of an array containing the parameters
param1.matrix = matrix(c(0.6,2,0.4,4,0.7,2,5,0.6,2,3,0.6,0.6,4,0.3,0.4,4), 
nrow = S, ncol = S, byrow = TRUE)
param2.matrix = matrix(c(0.8,0,0,2,0,5,0,0.8,0,0,0.8,0,4,0,0,4), nrow = S, 
ncol = S, byrow = TRUE)
param.array = array(c(param1.matrix, param2.matrix), c(S,S,2))
### Simulation of 3 sequences
seq3 = simulSM(E = alphabet, NbSeq = 3, lengthSeq = lengthSeq3, 
TypeSojournTime = "fij", init = vect.init, Ptrans = Pij, distr = distr.matrix, 
param = param.array, File.out = NULL)
                

################################
## Parametric simulation of several trajectories (3 trajectories of length 1000, 
## 10 000 and  2000 respectively),
## where the sojourn times depend only on the present state ("fi.")
## and the sojourn times are modelled by different distributions.
################################
## creation of the distribution matrix
distr.vect = c("dweibull", "pois", "geom", "nbinom")
## creation of an array containing the parameters
param.matrix = matrix(c(0.6,0.8,4,0,0.7,0,5,2), nrow = S, ncol = 2, byrow = TRUE)
### Simulation of 3 sequences without censoring
#seqFi3 = simulSM(E = alphabet, NbSeq = 3, lengthSeq = lengthSeq3, TypeSojournTime = "fi", 
#             init = vect.init, Ptrans = Pij, distr = distr.vect, param = param.matrix,
#             File.out = NULL)
### Simulation of 3 sequences with censoring at the beginning
#seqFi3 = simulSM(E = alphabet, NbSeq = 3, lengthSeq = lengthSeq3, TypeSojournTime = "fi", 
#             init = vect.init, Ptrans = Pij, distr = distr.vect, param = param.matrix, 
#             File.out = NULL, cens.beg = 1, cens.end = 0)
### Simulation of 3 sequences with censoring at the end
#seqFi3 = simulSM(E = alphabet, NbSeq = 3, lengthSeq = lengthSeq3, TypeSojournTime = "fi", 
#             init = vect.init, Ptrans = Pij, distr = distr.vect, param = param.matrix, 
#             File.out = NULL, cens.beg = 0, cens.end = 1)
### Simulation of 3 sequences with censoring at the beginning and at the end
#seqFi3 = simulSM(E = alphabet, NbSeq = 3, lengthSeq = lengthSeq3, TypeSojournTime = "fi", 
#             init = vect.init, Ptrans = Pij, distr = distr.vect, param = param.matrix,
#             File.out = NULL, cens.beg = 1, cens.end = 1)
                
################################
## Non-parametric simulation of several trajectories (3 trajectories of length 1000,  
## 10 000 and  2000 respectively),
## where the sojourn times depend on the present state and on the next state ("fij")
################################
## creation of a matrix corresponding to the conditionnal sojourn time distribution
Kmax = 4 
mat1 = matrix(c(0,0.5,0.4,0.6,0.3,0,0.5,0.4,0.7,0.2,0,0.3,0.4,0.6,0.2,0), nrow = S, 
ncol = S, byrow = TRUE)
mat2 = matrix(c(0,0.2,0.3,0.1,0.2,0,0.2,0.3,0.1,0.4,0,0.3,0.2,0.1,0.3,0), nrow = S, 
ncol = S, byrow = TRUE)
mat3 = matrix(c(0,0.1,0.3,0.1,0.3,0,0.1,0.2,0.1,0.2,0,0.3,0.3,0.3,0.4,0), nrow = S, 
ncol = S, byrow = TRUE)
mat4 = matrix(c(0,0.2,0,0.2,0.2,0,0.2,0.1,0.1,0.2,0,0.1,0.1,0,0.1,0), nrow = S, 
ncol = S, byrow = TRUE)
f <- array(c(mat1,mat2,mat3,mat4), c(S,S,Kmax))
### Simulation of 3 sequences
seqNP3 = simulSM(E = alphabet, NbSeq = 3, lengthSeq = lengthSeq3, TypeSojournTime = "fij", 
                init = vect.init, Ptrans = Pij, laws = f, File.out = NULL)


################################
## Non-parametric simulation of several trajectories (3 trajectories of length 1000, 
## 10 000 and  2000 respectively),
## where the sojourn times depend only on he next state ("fj")
################################
## creation of a matrix corresponding to the conditional sojourn time distributions
Kmax = 6
nparam.matrix = matrix(c(0.2,0.1,0.3,0.2,0.2,0,0.4,0.2,0.1,0,0.2,0.1,
 0.5,0.3,0.15,0.05,0,0,0.3,0.2,0.1,0.2,0.2,0), 
                       nrow = S, ncol = Kmax, byrow = TRUE)
### Simulation of 3 sequences without censoring
#seqNP3 = simulSM(E = alphabet, NbSeq = 3, lengthSeq = lengthSeq3, TypeSojournTime = "fj", 
#                init = vect.init, Ptrans = Pij, laws = nparam.matrix, File.out = NULL)
### Simulation of 3 sequences with censoring at the beginnig 
#seqNP3 = simulSM(E = alphabet, NbSeq = 3, lengthSeq = lengthSeq3, TypeSojournTime = "fj", 
#                init = vect.init, Ptrans = Pij, laws = nparam.matrix, File.out = NULL, 
#                cens.beg = 1, cens.end = 0)
### Simulation of 3 sequences with censoring at the end 
#seqNP3 = simulSM(E = alphabet, NbSeq = 3, lengthSeq = lengthSeq3, TypeSojournTime = "fj", 
#                init = vect.init, Ptrans = Pij, laws = nparam.matrix, File.out = NULL, 
#                cens.beg = 0, cens.end = 1)
### Simulation of 3 sequences with censoring at the beginning and at the end 
#seqNP3 = simulSM(E = alphabet, NbSeq = 3, lengthSeq = lengthSeq3, TypeSojournTime = "fj", 
#                init = vect.init, Ptrans = Pij, laws = nparam.matrix, File.out = NULL, 
#                cens.beg = 1, cens.end = 1)
}
}
\references{
VS Barbu, C Berard, D Cellier, M Sautreuil and N Vergne (2017), Parametric estimation of semi-Markov chains, submitted

}   
%R Varadhan and PD Gilbert (2009),  BB: An R Package for Solving a Large System of Nonlinear Equations and for Optimizing a High-Dimensional Nonlinear Objective Function, \emph{J. Statistical Software}, 32:4, \url{http://www.jstatsoft.org/v32/i04/}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Simulation}
\keyword{Semi-Markov models}% __ONLY ONE__ keyword per line
