simulMk <- function(E, nbSeq, lengthSeq, Ptrans, init, k, File.out = NULL){

  if( !is.vector(lengthSeq) ){
    stop("The parameter \"lenghtSeq\" must be a vector")
  }

  if( length(lengthSeq) != nbSeq ){
    stop("The length of \"lenghtSeq\" must be equal at the value of the parameter \"nbSeq\"")
  }

  if( !is.matrix(Ptrans) ){
    stop("The parameter \"Ptrans\" must be a matrix")
  }

  ifelse(rowSums(Ptrans) == 1, "", stop("The matrix \"Ptrans\" must be stochastic"))

  if ( sum(init) != 1 ){
    stop("The vector \"init\" must be equal to 1")
  }

  seq <- NULL
  y <- NULL
  for (n in 1:nbSeq) {
    for (i in 1:k){
      y[i]<- sample(E,1, prob=init) 
    }

    s<-length(E)
    for (i in 1:(lengthSeq[n]-k)){
      ind<- which(E==y[i+k-1])
      if (k > 1){
        for(j in (k-2):0) {ind<-ind+s^(j+1)*(which(E==y[i+j])-1) }  
      }
      y[i + k] <- sample(E, 1, prob = Ptrans[ind, ])
    }
    seq[[n]] = y

    ## Save in a fasta file
    if(!is.null(File.out)){
      if (file.exists(File.out)){
        write.fasta(sequences = y, names = "seq", file.out = File.out, open = "a")
      }else{
        write.fasta(y, names = "seq", file.out = File.out)
      }
    }
  }


  return(seq)
}
