\name{AIC_Mk}
\alias{AIC_Mk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
AIC (Markov model)
}
\description{
AIC
}
\usage{
AIC_Mk(seq, E, mu, Ptrans, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{seq}{
List of sequence(s)
}
  \item{E}{
Vector of state space 
}
  \item{mu}{
Vector of initial distribution
}
  \item{Ptrans}{
Matrix of transition probabilities
}
  \item{k}{
Order of Markov model
}
}
\details{
\eqn{AIC(M) = -2*log{L} + 2*M}, where L is the log-likelihood, M is the number of parameters of the model.
}
\value{
\item{AIC}{ List: value of AIC for each sequence}
}
\author{
Vlad Stefan Barbu, barbu@univ-rouen.fr \cr Caroline Berard, caroline.berard@univ-rouen.fr \cr Dominique Cellier, dominique.cellier@laposte.net \cr Mathilde Sautreuil, mathilde.sautreuil@etu.univ-rouen.fr \cr Nicolas Vergne, nicolas.vergne@univ-rouen.fr 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{simulSM},
\link{estimMk},
\link{simulMk},
\link{estimSM},
\link{LoglikelihoodSM},
\link{LoglikelihoodMk}
}
\examples{
alphabet = c("a","c","g","t")
S = length(alphabet)
# creation of the transition matrix
Pij = matrix(c(0,0.2,0.3,0.5,0.4,0,0.2,0.4,0.1,0.2,0,0.7,0.8,0.1,0.1,0), 
nrow = S, ncol = S, byrow = TRUE)

#Pij
#     [,1] [,2] [,3] [,4]
#[1,]  0.0  0.2  0.3  0.5
#[2,]  0.4  0.0  0.2  0.4
#[3,]  0.1  0.2  0.0  0.7
#[4,]  0.8  0.1  0.1  0.0


## Simulation of two sequences of length 20 and 50 respectively
seq2 = simulMk(E = alphabet, nbSeq = 2, lengthSeq = c(20,50), Ptrans = Pij, 
init = rep(1/4,4), k = 1)
                
#################################
## Computation of AIC
#################################
AIC_Mk(seq = seq2, E = alphabet, mu = rep(1/4,4), Ptrans = Pij, k = 1)

#[[1]]
#[1] 60.20263
#
#[[2]]
#[1] 115.7674
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{AIC}
\keyword{Markov models}% __ONLY ONE__ keyword per line
