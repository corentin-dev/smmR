InitialLawMk = function(E, seq, Ptrans, k){

  ## Verify the existence of the sequence 
  if( !is.list(seq) ){
    stop("The parameter \"seq\" should be a list")
  }

  ## get the state space size
  S = length(E)
  vect.seq<-NULL
  ## get the number of sequences
  nbSeq<-length(seq)

  if (nbSeq > 1 ){
    if (k == 1){
      init = .stationnary.law(Pest = Ptrans) ## !!! Warning if k > 1
      ## stationnary law
    }else{
      Nstart = as.vector(count(seq = unlist(seq), wordsize = 1, alphabet = E)) ## initial law k > 1 ???
      init = Nstart/sum(Nstart)
    }
  }else{
    if (k == 1){
      init = .stationnary.law(Pest = Ptrans) ## !!! Warning if k > 1
      ## stationnary law
    }else{
      Nstart = as.vector(count(seq = unlist(seq), wordsize = 1, alphabet = E)) ## initial law k > 1 ???
      init = Nstart/sum(Nstart)
    }
  }

  return( list(init = init) )
}
