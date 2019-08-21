### INPUTS
## seq: list of sequences
## alphabet: vector of state space 
## mu: vector of initial distribution
## puv : matrix of trasition probabilities 
## k : order of the Markov chain
LoglikelihoodMk = function(seq, E, mu, Ptrans, k){

  ## length of the state space
  S = length(E)
  if(dim(Ptrans)[1] != S || dim(Ptrans)[2] != S){
    stop("The size of the matrix Ptrans must be equal to SxS with S = length(E)")  
  }

  if( !is.matrix(Ptrans) ){
    stop("The parameter \"Ptrans\" must be a matrix")
  }
  
  ifelse(rowSums(Ptrans) == 1, "", stop("The matrix \"Ptrans\" must be stochastic"))
  
  if ( sum(mu) != 1 ){
    stop("The vector \"init\" must be equal to 1")
  }
  
  if(!is.list(seq)){
    stop("The parameter \"seq\" should be a list")
  }
  
  vect.seq<-NULL
  ## Get the number of sequences
  nbSeq<-length(seq)

  LNij = list()
  LNi = list()
  SeqP1 = NULL
  Ni = vector(length = S)
  ## Count the number of transitions from i to j in k steps
  Nij = matrix(0, nrow = S^k, ncol = S)
  for (i in 1:nbSeq){
    vect.seq = seq[[i]]
    Nij <- Nij + matrix(count(seq = vect.seq, wordsize = k+1, alphabet = E), byrow=TRUE, ncol=S)
    LNij[[i]] = matrix(count(seq = vect.seq, wordsize = k+1, alphabet = E), byrow=TRUE, ncol=S)
    ## Count the number of states i
    Ni <- Ni + as.vector(count(seq = vect.seq[1:(length(vect.seq)-k)], wordsize = k, alphabet = E))
    LNi[[i]] = as.vector(count(seq = vect.seq[1:(length(vect.seq)-k)], wordsize = k, alphabet = E))
    ## Get the first state
    SeqP1 = c(SeqP1, vect.seq[1])
  }
  
  lV = list()
  for (j in 1:nbSeq){
    s <- 0
    for (i in 1:k){ # Warning to initial law
      s <- s + log(mu[which(E==seq[[j]][i])])
    }  
    LNij.vect = LNij[[j]][which(Ptrans != 0)]
    Ptrans.vect = Ptrans[which(Ptrans != 0)]
    lV[[j]] <- s + sum(LNij.vect*log(Ptrans.vect))
  }
  
  # s<-as.numeric(s)     
  
  return (list(L = lV))
}
