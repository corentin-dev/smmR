BIC_Mk = function(seq, E, mu, Ptrans, k){
  
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
  ## Kpar: number of parameters of the model
  Kpar = (S-1)*S^k
  
  ## Computation of the log-likelihood for all sequences
  res = LoglikelihoodMk(seq = seq, E = E, mu = mu, Ptrans = Ptrans, k = k)
  LV = res$L
  nbSeq = length(seq)
  
  BIC.list = list()
  for (k in 1:nbSeq) {
    n = length(seq[[k]])
    BIC.list[[k]] = -2*LV[[k]] + log(n)*Kpar 
  } 
  
  
  return(BIC = BIC.list)
}
