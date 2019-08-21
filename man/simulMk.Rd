\name{simulMk}
\alias{simulMk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of a k-th order Markov chain
}
\description{
Simulation of a k-th order Markov chain starting from a transition matrix and an initial distribution.
}
\usage{
simulMk(E, nbSeq, lengthSeq, Ptrans, init, k, File.out = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{
Vector of state space
}
  \item{nbSeq}{
Number of simulated sequences 
}
  \item{lengthSeq}{
Vector of size nbSeq containing the lengths of each simulated sequence
}
  \item{Ptrans}{
Matrix of transition probabilities of size (S^(k))xS, with S = length(E) 
}
  \item{init}{
Vector of initial distribution of length S
}
  \item{k}{
Order of the Markov chain
}
  \item{File.out}{
Name of the fasta file for saving the sequences.  
If File.out = NULL, no file is created
}
}
\details{
The sizes of init and Ptrans depend on S, the length of E.
The rows of the transition matrix sums to 1.

For k=1, the transition matrix is defined by \eqn{Ptrans(i,j) = P(X_{m+1} = j | X_m = i)} and 
the initial distribution is \eqn{init = P(X_1 = i)}. For k > 1 we have similar expressions.


The first element of lengthSeq corresponds to the length of the first sequence and so on.
}
\value{
simulMk returns a list of sequences of size lengthSeq simulated with a k-th order Markov chain of parameters init and Ptrans with state space E.

If the parameter File.out is not equal to NULL, a file in fasta format containing the sequence(s) will be created.
}

%%\references{
%%}


\author{
Vlad Stefan Barbu, barbu@univ-rouen.fr \cr Caroline Berard, caroline.berard@univ-rouen.fr \cr Dominique Cellier, dominique.cellier@laposte.net \cr Mathilde Sautreuil, mathilde.sautreuil@etu.univ-rouen.fr \cr Nicolas Vergne, nicolas.vergne@univ-rouen.fr 
}


%%\note{
%%
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
estimMk,
simulSM,
estimSM
}

\examples{

### Example 1 ###
# Second order model with state space {a,c,g,t}
E <- c("a","c","g","t")
S = length(E)
init.distribution <- c(1/6,1/6,5/12,3/12)
k<-2
p <- matrix(0.25, nrow = S^k, ncol = S)

# We simulate 3 sequences of size 1000, 10000 and 2000 respectively.
simulMk(E = E, nbSeq = 3, lengthSeq = c(1000, 10000, 2000), Ptrans = p, 
init = init.distribution, k = k)

### Example 2 ###
# first order model with state space {1,2,3}
E <- c(1,2,3)
S <- length(E)
init.distr <- rep(1/S, 3)
p <- matrix(c(0.3,0.2,0.5,0.1,0.6,0.3,0.2,0.4,0.4), nrow = 3, byrow = TRUE)

# We simulate one sequence of size 100
simulMk(E = E, nbSeq = 1, lengthSeq = 100, Ptrans = p, init = init.distr, k = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Markov models}
\keyword{Simulation}% __ONLY ONE__ keyword per line
