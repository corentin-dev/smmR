\name{estimSM}
\alias{estimSM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of a semi-Markov chain
}
\description{
Estimation of a semi-Markov chain starting from one or several sequences. This estimation can be parametric or non-parametric, non-censored, censored at the beginning and/or at the end of the sequence, with one or several trajectories. Several parametric distributions are considered (Uniform, Geometric, Poisson, Discrete Weibull and Negative Binomial). 
}
\usage{
# parametric case
estimSM(file = NULL, seq, E, TypeSojournTime = "fij", distr, cens.end = 0, 
cens.beg = 0)

estimSM(file, seq = NULL, E, TypeSojournTime = "fij", distr, cens.end = 0, 
cens.beg = 0)

# non-parametric case
estimSM(file = NULL, seq, E, TypeSojournTime = "fij", distr = "NP", cens.end = 0, 
cens.beg = 0)

estimSM(file, seq = NULL, E, TypeSojournTime = "fij", distr = "NP", cens.end = 0, 
cens.beg = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
Path of the file in fasta format which contains the sequences from which to estimate  
}
  \item{seq}{
List of the sequence(s)
}
  \item{E}{
Vector of state space of length S
}
  \item{TypeSojournTime}{
Character: "fij", "fi", "fj", "f" (for more explanations, see Details)
}
  \item{distr}{
- "NP" for nonparametric case, laws have to be used, param is useless

- Matrix of distributions of size SxS if TypeSojournTime is equal to "fij";

- Vector of distributions of size S if TypeSojournTime is equal to "fi" or "fj";

- A distribution if TypeSojournTime is equal to "f".

The distributions to be used in distr must be one of "uniform", "geom", "pois", "dweibull", "nbinom".
}
  \item{cens.beg}{
Type of censoring at the beginning of sample paths;  1 (if the first sojourn time is censored) or 0 (if not). All the sequences must be censored in
the same way.}
  \item{cens.end}{
Type of censoring at the end of sample paths;  1 (if the last sojourn time is censored) or 0 (if not). All the sequences must be censored in
the same way.
}
}
\details{
This function estimates a semi-Markov model in parametric and non-parametric case, taking into account the type of sojourn time and the censoring described in references.
The non-parametric estimation concerns sojourn time distributions defined by the user.
For the parametric estimation, several discrete distributions are considered (see below).

The difference between the Markov model and the semi-Markov model concerns the modelisation of the sojourn time. With a Markov chain, the sojourn time distribution is modeled by a Geometric distribution (in discrete time). With a semi-Markov chain, the sojourn time can be any arbitrary distribution.
In this package, the available distribution for a semi-Markov model are :
\itemize{
\item Uniform: \eqn{f(x) = 1/n} for \eqn{a \le x \le b}, with \eqn{n = b-a+1}
\item Geometric: \eqn{f(x) = \theta (1-\theta)^x} for \eqn{x = 0, 1, 2,\ldots,n}, \eqn{0 < \theta \le 1}, with \eqn{n > 0} and \eqn{\theta} is the probability of success
\item Poisson: \eqn{f(x) = (\lambda^x exp(-\lambda))/x!} for \eqn{x = 0, 1, 2,\ldots,n}, with \eqn{n > 0} and \eqn{\lambda > 0}
\item Discrete Weibull of type 1: \eqn{f(x)=q^{(x-1)^{\beta}}-q^{x^{\beta}}, x=1,2,3,\ldots,n}, with \eqn{n > 0}, \eqn{q} is the first parameter and \eqn{\beta} is the second parameter
\item Negative binomial: \eqn{f(x)=\Gamma(x+\alpha)/(\Gamma(\alpha) x!) (\alpha/(\alpha+\mu))^{\alpha} (\mu/(\alpha+\mu))^x}, for \eqn{x = 0, 1, 2,\ldots,n}, \eqn{n > 0,\ \Gamma } is the Gamma function, \eqn{\alpha} is the parameter of overdispersion and \eqn{\mu} is the mean
\item Non-parametric
}
We define :
\itemize{
\item the semi-Markov kernel \eqn{q_{ij}(k) = P( J_{m+1} = j, T_{m+1} - T_{m} = k | J_{m} = i )};
\item the transition matrix \eqn{(Ptrans(i,j))_{i,j} \in E} of the embedded Markov chain \eqn{J = (J_m)_m}, \eqn{Ptrans(i,j) = P( J_{m+1} = j | J_m = i )};
\item the initial distribution \eqn{init = P(J_1 = i) = P(Y_1 = i)};
\item the conditional sojourn time distributions \eqn{(f_{ij}(k))_{i,j} \in E,\ k \in N ,\ f_{ij}(k) = P(T_{m+1} - T_m = k | J_m = i, J_{m+1} = j )}, f is specified by the argument "param" in the parametric case and by "laws" in the non-parametric case.
}

The estimation of the transition matrix of the embedded Markov chain is \eqn{\widehat{Ptrans}(i,j) = N_{ij}/N_{i.}.}

The estimation of the initial law is the limit law if the number of sequences is equal to 1, else the estimation of the inital law is 
\eqn{init = N_i^l/N^l}, where \eqn{N_i^l} is the number of times the state i appears in all the sequences and \eqn{N^l} is the size of sequences.

In the parametric case, the distribution of sojourn time is calculated with the estimated parameters.

Note that \eqn{q_{ij}(k) = Ptrans(i,j)*f_{ij}(k)} in the general case (depending on the present state and on the next state). For particular cases, we replace \eqn{f_{ij}(k)} by \eqn{f_{i.}(k)} (depending on the present state "fi."), \eqn{f_{.j}(k)} (depending on the next state "f.j") and \eqn{f_{..}(k)} (depending neither on the present state nor on the next state "f").

In this package we can choose differents types of sojourn time. Four options are available for the sojourn times:\itemize{
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

At any use of the function, a file in .txt format is created. This file contains the output of the function.

}
\value{
\item{init}{Vector of initial distribution}
\item{Ptrans}{Matrix of transition probabilities of the embedded Markov chain \eqn{J=(J_m)_{m}} of size SxS}
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
\item{q}{Array of size SxSxKmax with the estimated semi-Markov kernel}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}


\author{
Vlad Stefan Barbu, barbu@univ-rouen.fr \cr Caroline Berard, caroline.berard@univ-rouen.fr \cr Dominique Cellier, dominique.cellier@laposte.net \cr Mathilde Sautreuil, mathilde.sautreuil@etu.univ-rouen.fr \cr Nicolas Vergne, nicolas.vergne@univ-rouen.fr 
}


%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{simulSM},
\link{estimMk},
\link{simulMk}
}
\examples{
alphabet = c("a","c","g","t")
S = length(alphabet)
# creation of the transition matrix
Pij = matrix(c(0,0.2,0.3,0.5,0.4,0,0.2,0.4,0.1,0.2,0,0.7,0.8,0.1,0.1,0), nrow = S, 
ncol = S, byrow = TRUE)

Pij
#     [,1] [,2] [,3] [,4]
#[1,]  0.0  0.2  0.3  0.5
#[2,]  0.4  0.0  0.2  0.4
#[3,]  0.1  0.2  0.0  0.7
#[4,]  0.8  0.1  0.1  0.0

################################
## Parametric estimation of a trajectory (of length equal to 5000)
## where the sojourn times depend neither on the present state nor on the next state ("f")
################################
## Simulation of a sequence of length 5000
seq5000 = simulSM(E = alphabet, NbSeq = 1, lengthSeq = 5000, TypeSojournTime = "f",
               init = c(1/4,1/4,1/4,1/4), Ptrans = Pij, distr = "pois", param = 2)


## Estimation of the simulated sequence
estSeq5000 = estimSM(seq = seq5000, E = alphabet, TypeSojournTime = "f",
                  distr = "pois", cens.end = 0, cens.beg = 0)

# initial distribution estimated
estSeq5000$init
# [1] 0.3592058 0.1456077 0.1600481 0.3351384

# transition matrix estimated
estSeq5000$Ptrans
#          [,1]       [,2]       [,3]      [,4]
#[1,] 0.0000000 0.21775544 0.30150754 0.4807370
#[2,] 0.4297521 0.00000000 0.18181818 0.3884298
#[3,] 0.1052632 0.23308271 0.00000000 0.6616541
#[4,] 0.8348294 0.08976661 0.07540395 0.0000000

# estimated parameter
estSeq5000$param
# [1] 2.007822 0.000000

# estimated semi-Markov kernel
estSeq5000$q
# , , 1
#           [,1]       [,2]       [,3]       [,4]
#[1,] 0.00000000 0.02924038 0.04048668 0.06455377
#[2,] 0.05770747 0.00000000 0.02441470 0.05215867
#[3,] 0.01413482 0.03129854 0.00000000 0.08884747
#[4,] 0.11210159 0.01205393 0.01012531 0.00000000
#, , 2
#           [,1]       [,2]       [,3]      [,4]
#[1,] 0.00000000 0.05870948 0.08129005 0.1296125
#[2,] 0.11586631 0.00000000 0.04902036 0.1047253
#[3,] 0.02838021 0.06284189 0.00000000 0.1783899
#[4,] 0.22508003 0.02420215 0.02032981 0.0000000
#, , 3
#          [,1]       [,2]       [,3]      [,4]
#[1,] 0.0000000 0.05893909 0.08160797 0.1301194
#[2,] 0.1163195 0.00000000 0.04921208 0.1051349
#[3,] 0.0284912 0.06308767 0.00000000 0.1790876
#[4,] 0.2259603 0.02429681 0.02040932 0.0000000
#, , 4
#           [,1]       [,2]       [,3]       [,4]
#[1,] 0.00000000 0.03944640 0.05461809 0.08708551
#[2,] 0.07784959 0.00000000 0.03293636 0.07036405
#[3,] 0.01906842 0.04222293 0.00000000 0.11985865
#[4,] 0.15122935 0.01626122 0.01365943 0.00000000
#, , 5
#            [,1]        [,2]        [,3]       [,4]
#[1,] 0.000000000 0.019800336 0.027415850 0.04371305
#[2,] 0.039077026 0.000000000 0.016532588 0.03531962
#[3,] 0.009571498 0.021194032 0.000000000 0.06016370
#[4,] 0.075910402 0.008162409 0.006856423 0.00000000
#, , 6
#            [,1]        [,2]        [,3]      [,4]
#[1,] 0.000000000 0.007951110 0.011009229 0.0175536
#[2,] 0.015691942 0.000000000 0.006638898 0.0141831
#[3,] 0.003843573 0.008510768 0.000000000 0.0241596
#[4,] 0.030482913 0.003277733 0.002753295 0.0000000
#, , 7
#            [,1]        [,2]         [,3]        [,4]
#[1,] 0.000000000 0.002660735 0.0036840950 0.005874085
#[2,] 0.005251104 0.000000000 0.0022216210 0.004746190
#[3,] 0.001286202 0.002848018 0.0000000000 0.008084696
#[4,] 0.010200710 0.001096851 0.0009213545 0.000000000
#, , 8
#             [,1]         [,2]         [,3]        [,4]
#[1,] 0.0000000000 0.0007631832 0.0010567152 0.001684874
#[2,] 0.0015061831 0.0000000000 0.0006372313 0.001361358
#[3,] 0.0003689234 0.0008169018 0.0000000000 0.002318947
#[4,] 0.0029258870 0.0003146115 0.0002642737 0.000000000
#, , 9
#             [,1]         [,2]         [,3]         [,4]
#[1,] 0.000000e+00 1.915420e-04 2.652120e-04 0.0004228658
#[2,] 3.780184e-04 0.000000e+00 1.599309e-04 0.0003416705
#[3,] 9.259156e-05 2.050242e-04 0.000000e+00 0.0005820041
#[4,] 7.343325e-04 7.896048e-05 6.632681e-05 0.0000000000

\donttest{
################################
## Parametric estimation of a trajectory (of length equal to 5000), 
## where sojourn times do not depend neither on the present state nor on the next state
## and the sequence is censored at the beginning.
################################
## Simulation of a sequence of length 5000
#seq5000 = simulSM(E = alphabet, NbSeq = 1, lengthSeq = 5000, TypeSojournTime = "f", 
#                init = c(1/4,1/4,1/4,1/4), Ptrans = Pij, distr = "pois", param = 2, 
#                cens.beg = 1, cens.end = 0)
                
                
## Estimation of the simulated sequence                
#estSeq5000 = estimSM(seq = seq5000, E = alphabet, TypeSojournTime = "f", 
#                   distr = "pois", cens.end = 0, cens.beg = 1)


################################
## Parametric estimation of a trajectory (of length equal to 5000), 
## where sojourn times do not depend neither on the present state nor on the next state
## and the sequence is censored at the beginning and at the end
################################
## Simulation of a sequence of legnth 5000
#seq5000 = simulSM(E = alphabet, NbSeq = 1, lengthSeq = 5000, TypeSojournTime = "f", 
#                init = c(1/4,1/4,1/4,1/4), Ptrans = Pij, distr = "pois", param = 2, 
#                cens.beg = 1, cens.end = 1)
                
                
## Estimation of the simulated sequence                
#estSeq5000 = estimSM(seq = seq5000, E = alphabet, TypeSojournTime = "f", 
#                   distr = "pois", cens.end = 1, cens.beg = 1)
                   
################################
## Parametric simulation of several trajectories (3 trajectories of length 1000, 10 000 
## and 2000 respectively),
## where the sojourn times depend on the present state and on the next state
## and the sojourn time distributions are modeled by different distributions.
################################
lengthSeq3 = c(1000, 10000, 2000)
## creation of the initial distribution
vect.init = c(1/4,1/4,1/4,1/4)
## creation of the distribution matrix
distr.matrix = matrix(c("dweibull", "pois", "geom", "nbinom", "geom", "nbinom", 
"pois", "dweibull", "pois", "pois", "dweibull", "geom", "pois","geom", "geom", 
"nbinom"), nrow = S, ncol = S, byrow = TRUE)
## creation of an array containing the parameters
param1.matrix = matrix(c(0.6,2,0.4,4,0.7,2,5,0.6,2,3,0.6,0.6,4,0.3,0.4,4), nrow = S, 
ncol = S, byrow = TRUE)
param2.matrix = matrix(c(0.8,0,0,2,0,5,0,0.8,0,0,0.8,0,4,0,0,4), nrow = S, ncol = S, 
byrow = TRUE)
param.array = array(c(param1.matrix, param2.matrix), c(S,S,2))
### Simulation of 3 sequences
seq3 = simulSM(E = alphabet, NbSeq = 3, lengthSeq = lengthSeq3, TypeSojournTime = "fij", 
init = vect.init, Ptrans = Pij, distr = distr.matrix, param = param.array, File.out = NULL)
                
## Estimation of the simulated sequence                
estSeq3 = estimSM(seq = seq3, E = alphabet, TypeSojournTime = "fij", 
                   distr = distr.matrix, cens.end = 0, cens.beg = 0)
                   
################################
## Non-Parametric simulation of several trajectories (3 trajectories of length 1000, 10 000 
## and 2000 respectively),
## where the sojourn times depend on the present state and on the next state
################################
lengthSeq3 = c(1000, 10000, 2000)
## creation of the initial distribution
vect.init = c(1/4,1/4,1/4,1/4)
## creation of an array containing the conditional distributions
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
seq3.NP = simulSM(E = alphabet, NbSeq = 3, lengthSeq = lengthSeq3, 
TypeSojournTime = "fij", init = vect.init, Ptrans = Pij, distr = "NP", laws = f, 
File.out = NULL)
                
## Estimation of the simulated sequence                
estSeq3.NP = estimSM(seq = seq3.NP, E = alphabet, TypeSojournTime = "fij", 
                   distr = "NP", cens.end = 0, cens.beg = 0)

                 
# initial distribution estimated
estSeq3.NP$init
# [1] 0.1856190 0.2409524 0.2975714 0.2758571

# transition matrix estimated
estSeq3.NP$Ptrans
#            [,1]       [,2]       [,3]      [,4]
# [1,] 0.00000000 0.20560325 0.29191143 0.5024853
# [2,] 0.40614334 0.00000000 0.19795222 0.3959044
# [3,] 0.08932039 0.21941748 0.00000000 0.6912621
# [4,] 0.81206817 0.09120221 0.09672962 0.0000000

# parameter estimated
estSeq3.NP$laws
#, , 1
#          [,1]      [,2]      [,3]      [,4]
# [1,] 0.0000000 0.4769231 0.4009288 0.6016187
# [2,] 0.3053221 0.0000000 0.4540230 0.3764368
# [3,] 0.7717391 0.2389381 0.0000000 0.2794944
# [4,] 0.4140669 0.5959596 0.2000000 0.0000000
# , , 2
#           [,1]      [,2]      [,3]      [,4]
# [1,] 0.00000000 0.2175824 0.2801858 0.1052158
# [2,] 0.16526611 0.0000000 0.2356322 0.2959770
# [3,] 0.07608696 0.3716814 0.0000000 0.3089888
# [4,] 0.18774816 0.1161616 0.3142857 0.0000000
# , , 3
#           [,1]       [,2]      [,3]       [,4]
# [1,] 0.00000000 0.09450549 0.3188854 0.08363309
# [2,] 0.31092437 0.00000000 0.1034483 0.19540230
# [3,] 0.06521739 0.20353982 0.0000000 0.32022472
# [4,] 0.29892229 0.28787879 0.3666667 0.00000000
#, , 4
#           [,1]      [,2]      [,3]       [,4]
# [1,] 0.00000000 0.2109890 0.0000000 0.20953237
# [2,] 0.21848739 0.0000000 0.2068966 0.13218391
# [3,] 0.08695652 0.1858407 0.0000000 0.09129213
# [4,] 0.09926262 0.0000000 0.1190476 0.00000000

# semi-Markovian kernel estimated
estSeq3.NP$q
# , , 1
#           [,1]       [,2]       [,3]      [,4]
# [1,] 0.00000000 0.09805694 0.11703570 0.3023046
# [2,] 0.12400455 0.00000000 0.08987486 0.1490330
# [3,] 0.06893204 0.05242718 0.00000000 0.1932039
# [4,] 0.33625058 0.05435283 0.01934592 0.0000000
# , , 2
#            [,1]       [,2]       [,3]       [,4]
# [1,] 0.000000000 0.04473565 0.08178943 0.05286941
# [2,] 0.067121729 0.00000000 0.04664391 0.11717861
# [3,] 0.006796117 0.08155340 0.00000000 0.21359223
# [4,] 0.152464302 0.01059420 0.03040074 0.00000000
#, , 3
#            [,1]       [,2]       [,3]       [,4]
# [1,] 0.000000000 0.01943064 0.09308631 0.04202440
# [2,] 0.126279863 0.00000000 0.02047782 0.07736064
# [3,] 0.005825243 0.04466019 0.00000000 0.22135922
# [4,] 0.242745279 0.02625518 0.03546753 0.00000000
#, , 4
#           [,1]       [,2]       [,3]      [,4]
# [1,] 0.00000000 0.04338003 0.00000000 0.1052869
# [2,] 0.08873720 0.00000000 0.04095563 0.0523322
# [3,] 0.00776699 0.04077670 0.00000000 0.0631068
# [4,] 0.08060801 0.00000000 0.01151543 0.0000000

#---------------------------------------------#
alphabet = c("0","1")
S = length(alphabet)
# creation of the transition matrix
Pij = matrix(c(0,1,1,0), nrow = S, ncol = S, byrow = TRUE)
distr = matrix(c("nbinom", "pois", "geom", "geom"), nrow = S, ncol = S, byrow = TRUE)
param = array(c(matrix(c(2,5,0.4,0.7), nrow = S, ncol = S, byrow = TRUE), matrix(c(6,0,0,0), 
nrow = S, ncol = S, byrow = TRUE)), c(S,S,2)) 

################################
## Parametric estimation of a trajectory (of length  equal to 5000) 
## where the state space is {"0","1"}
################################
## Simulation of a sequence of length 5000
seq2 = simulSM(E = alphabet, NbSeq = 2, lengthSeq = c(5000,1000), TypeSojournTime = "fij", 
                init = c(1/2,1/2), Ptrans = Pij, distr = distr, param = param)
                
                
## Estimation of the simulated sequence                
estSeq2 = estimSM(seq = seq2, E = alphabet, TypeSojournTime = "fij", 
                   distr = distr, cens.end = 1, cens.beg = 1)


#---------------------------------------------#

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Estimation}
\keyword{Semi-Markov models}% __ONLY ONE__ keyword per line
