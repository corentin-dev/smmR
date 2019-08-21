\name{BIC_Mk}
\alias{BIC_Mk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BIC (Markov model)
}
\description{
BIC
}
\usage{
BIC_Mk(seq, E, mu, Ptrans, k)
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
Order of the Markov chain
}
}
\details{
\eqn{BIC(M) = -2*log{L} + log(n)*M}, where L is the log-likelihood, M is the number of parameters of the model and n is the size of the sequence.
}
\value{
\item{BIC}{ List: value of BIC for each sequence}
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
\link{LoglikelihoodMk},
\link{AIC_Mk}
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
seq2 = simulMk(E = alphabet, nbSeq = 2, lengthSeq = c(20,50), 
Ptrans = Pij, init = rep(1/4,4), k = 1)
                
#################################
## Computation of BIC
#################################
BIC_Mk(seq = seq2, E = alphabet, mu = rep(1/4,4), Ptrans = Pij, k = 1)

#[[1]]
#[1] 78.39401
#
#[[2]]
#[1] 133.7015
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BIC}
\keyword{Markov models}% __ONLY ONE__ keyword per lin
