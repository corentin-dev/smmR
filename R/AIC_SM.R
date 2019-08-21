AIC_SM = function(seq, E, mu, Ptrans, distr = "NP", param = NULL, laws = NULL, TypeSojournTime){
  
  S = length(E)
  if(dim(Ptrans)[1] != S || dim(Ptrans)[2] != S){
    stop("The size of the matrix Ptrans must be equal to SxS with S = length(E)")  
  }
  
  if( (distr == "NP" && is.null(laws)) || (length(distr) == 1 && distr != "NP" && !is.null(laws)) ){
    stop("The parameter \"param\" must be used with the parameter \"distr\"")
  }
  
  if ( (is.matrix(distr) || is.array(distr) || distr != "NP") && is.null(param) ){
    stop("The parameter \"param\" must be used with the parameter \"distr\"")
  }
  
  if( is.null(param) && is.null(laws) ){
    stop("One of the two parameters \"param\" (with distr) or \"laws\" (without distr) should be given")
  }
  
  if(!is.list(seq)){
    stop("The parameter \"seq\" should be a list")
  }
  
  ## All sequences 
  nbSeq<-length(seq)

  ## TypeSojournTime
  TYPES<-c("fij", "fi", "fj", "f")
  type<-pmatch(TypeSojournTime,TYPES)
  
  Compte = 0 
  if ( type == 1 ){
    
    if( length(distr) == 1 && distr == "NP" ){
      Kmax = dim(laws)[3]
      Kpar = (Kmax-1)*S*(S-2)
    }else{
      diago = seq(from = 1, to = S*S, by = S+1)
      lois = distr[-diago]
      
      g = which(lois == "geom")
      Ng = length(g)
      
      u = which(lois == "unif")
      Nu = length(u)
      
      p = which(lois == "pois")
      Np = length(p)
      
      nb = which(lois == "nbinom")
      Nnb = length(nb)
      
      dw = which(lois == "dweibull")
      Ndw = length(dw)
      
      Kpar = Ng + Np + Nu + 2*Ndw + 2*Nnb
    }
    
  }else if ( type == 2 ){
    
    if( length(distr) == 1 && distr == "NP" ){
      Kmax = dim(laws)[3]
      Kpar = (Kmax-1)*S*(S-2)
    }else{
      lois = distr
      
      g = which(lois == "geom")
      Ng = length(g)
      
      u = which(lois == "unif")
      Nu = length(u)
      
      p = which(lois == "pois")
      Np = length(p)
      
      nb = which(lois == "nbinom")
      Nnb = length(nb)
      
      dw = which(lois == "dweibull")
      Ndw = length(dw)
      
      Kpar = Ng + Np + Nu + 2*Ndw + 2*Nnb
    }
    
  }else if ( type == 3 ){
    
    if( length(distr) == 1 && distr == "NP" ){
      Kmax = dim(laws)[3]
      Kpar = (Kmax-1)*S*(S-2)
    }else{
      lois = distr
      
      g = which(lois == "geom")
      Ng = length(g)
      
      u = which(lois == "unif")
      Nu = length(u)
      
      p = which(lois == "pois")
      Np = length(p)
      
      nb = which(lois == "nbinom")
      Nnb = length(nb)
      
      dw = which(lois == "dweibull")
      Ndw = length(dw)
      
      Kpar = Ng + Np + Nu + 2*Ndw + 2*Nnb
    }    
  
  }else if ( type == 4 ){
    
    if( length(distr) == 1 && distr == "NP" ){
      Kmax = dim(laws)[3]
      Kpar = (Kmax-1)*S*(S-2)
    }else{
      lois = distr
      
      g = which(lois == "geom")
      Ng = length(g)
      
      u = which(lois == "unif")
      Nu = length(u)
      
      p = which(lois == "pois")
      Np = length(p)
      
      nb = which(lois == "nbinom")
      Nnb = length(nb)
      
      dw = which(lois == "dweibull")
      Ndw = length(dw)
      
      Kpar = Ng + Np + Nu + 2*Ndw + 2*Nnb
    }    
  
  }else{
    stop("The parameter \"TypeSojournTime\" must be one of \"fij\", \"fi\", \"fj\" et \"f\".")    
  }
  
  ## Computation of the log-likelihood for all  sequences 
  res = LoglikelihoodSM(seq = seq, E = E, mu = mu, Ptrans = Ptrans, distr = distr, param = param, 
                        laws = laws, TypeSojournTime = TypeSojournTime)
  LV = res$L
  
  AIC.list = list()
  for (k in 1:nbSeq) {
    AIC.list[[k]] = -2*LV[[k]] + 2*Kpar 
  } 
  
  
  return(AIC = AIC.list)
  
}
