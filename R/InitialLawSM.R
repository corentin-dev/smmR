InitialLawSM = function(E, seq, q){
  
  Ptrans = apply(q, c(1,2), sum)
  
  if(!is.list(seq)){
    stop("the parameter seq should be a list")
  }
  
  ## Size of the state space
  S<-length(E)
  
  ## All sequences 
  nbSeq<-length(seq)
  Kmax = 0
  J<-list()
  T<-list()
  L<-list()
  M<-list()
  
  for (k in 1:nbSeq){
    ## Manipulations
    JT<-.donneJT(unlist(seq[[k]]))
    J[[k]]<-.code(JT[1,],E)
    T[[k]]<-as.numeric(JT[2,])
    L[[k]]<-T[[k]] - c(1,T[[k]][-length(T[[k]]-1)]) ## sojourn times
    Kmax <- max(Kmax, max(L[[k]])) 
  }
  
  
  ## J: list
  ## L: list
  ## S: state space size
  ## Kmax: maximal sojourn time
  ################
  ## get the counts
  res<-.comptage(J,L,S,Kmax)
  
  Nstart = res$Nstarti
  #########
  
  if(nbSeq > 1){
    init1 = Nstart/sum(Nstart)
  }else{
    ## Computation of the limit law
    init1 = .limit.law(q = q, pij = Ptrans)
  }
  
  return( list(init = init1) )
}
