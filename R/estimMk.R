estimMk<-function(file = NULL, seq = NULL, E, k){  

  #require(seqinr)

  ## get the sequences from the given file
  if(!is.null(file)){
    fasta = read.fasta(file = file)
    seq = getSequence(fasta)
    attr = getAnnot(fasta)
  }

  ## verify the existence of the sequence
  if( !is.list(seq) ){
    stop("The parameter \"seq\" should be a list")
  }

  ## get the state space size
  S = length(E)
  vect.seq<-NULL
  ## get the number of sequences
  nbSeq<-length(seq)

  SeqP1 = NULL
  Ni = vector(length = S)
  ## count the number of transitions from i to j in k steps
  Nij = matrix(0, nrow = S^k, ncol = S)
  for (i in 1:nbSeq){
    vect.seq = seq[[i]]
    Nij <- Nij + matrix(count(seq = vect.seq, wordsize = k+1, alphabet = E), byrow=TRUE, ncol=S)
    ## count the number of states i
    Ni <- Ni + as.vector(count(seq = vect.seq[1:(length(vect.seq)-k)], wordsize = k, alphabet = E))
    ## get the state at position 1 for each sequence
    SeqP1 = c(SeqP1, vect.seq[1])
  }


  ## Verify the existence of all transitions
  if( length(which(Nij == 0)) > 0 ){
    warning("missing transitions")
  }

  ## Verify the existence of all states
  if( length(which(Ni == 0)) > 0 ){
    warning("missing observed states")
  }

  ## Compute the transition matrix
  Ptrans<-Nij/(Ni)
  Ptrans[which(Ni == 0),] = rep(0, S)

  if (nbSeq > 1 ){
    if (k == 1){
      init = .stationnary.law(Pest = Ptrans) ## !!!Warning if k > 1
      ## stationnary law
    }else{
      Nstart = as.vector(count(seq = unlist(seq), wordsize = 1, alphabet = E)) ## initial law for k > 1 ???
      init = Nstart/sum(Nstart)
    }
  }else{
    if (k == 1){
      init = .stationnary.law(Pest = Ptrans) ## !!! Warning if k > 1
      ## stationnary law
    }else{
      Nstart = as.vector(count(seq = unlist(seq), wordsize = 1, alphabet = E)) ## initial law for k > 1 ???
      init = Nstart/sum(Nstart)
    }
  }
  return(list(init = init, Ptrans = Ptrans))
}
