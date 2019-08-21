####################
#### INPUT :
## E: alphabet
## lengthSeq: lengths of sequences
## TypeSojournTime: "fij", "fi", "fj", "f"
## init: initial law
## Ptrans: transition matrix
## distr: matrix of distributions or vector of distributions or a distribution according to the type of sojourn times
## param: - array (in param[,,1] the first parameter of the distribution and in param[,,2] the second) if fij
##        - matrix (in param[,1] the first parameter of the distribution and in param[,2] the second) if fi ou fj
##        - vector (in param[1] the first parameter of the distribution and in param[2] the second) if f
## cens.beg : - 1 if censoring at the beginning
##            - 0 if not
## cens.end : - 1 if censoring at the end
##            - 0 if not
####################
#### OUTPUT : simulated sequence
simulSM<-function(E, NbSeq, lengthSeq, TypeSojournTime = "fij", init, Ptrans, distr = "NP", param = NULL, laws = NULL, cens.beg = 0, cens.end = 0, File.out = NULL){
  
  #require(DiscreteWeibull)
  
  S = length(E)
  if(dim(Ptrans)[1] != S || dim(Ptrans)[2] != S){
    stop("The size of the matrix Ptrans must be equal to SxS with S = length(E)")  
  }
  
  if( length(lengthSeq) != NbSeq ){
    stop("The length of \"lengthSeq\" must be equal at the value of NbSeq")  
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

  ## TypeSojournTime
  TYPES<-c("fij", "fi", "fj", "f")
  type<-pmatch(TypeSojournTime,TYPES)
  
  out = list()
  for (m in 1:NbSeq){
    J<-NULL
    J[1]<- sample(E,1, prob=init)
    T<-NULL
    i<-1
    s = 1
    t<-1
    if(length(distr) == 1 && distr == "NP"){
      
      f<-laws ## non-parametric distributions
      while (t <= lengthSeq[m]) {
        J[i + 1] <- sample(E, 1, prob = Ptrans[which(E==J[i]), ])
        ## fij
        if(type == 1){
          
          if( !is.array(f) ){
            stop("laws must be an array")
          }
          
          if( dim(f)[1] != S || dim(f)[2] != S ){
            stop("laws must be an array of size SxSxKmax")            
          }
          
          Kmax<-dim(f)[3]
          ## change code ?
          T[i] <- t + sample(1:Kmax,1,prob = f[.code(J[i],E), .code(J[i+1],E),])
          t<-T[i]
          i<- i+1
          
          ## fi  
        }else if(type == 2){
          
          if( dim(f)[1] != S ){
            stop("laws must be a matrix of size SxKmax")            
          }
          
          Kmax<-dim(f)[2]
          ## change code ?
          T[i] <- t + sample(1:Kmax,1,prob = f[.code(J[i],E),])
          t<-T[i]
          i<- i+1
          
          ## fj  
        }else if(type == 3){
          
          if( dim(f)[1] != S ){
            stop("laws must be a matrix of size SxKmax")            
          }
          
          Kmax<-dim(f)[2]
          T[i] <- t + sample(1:Kmax,1,prob = f[.code(J[i+1],E),])
          t<-T[i]
          i<- i+1
          
          ## f
        }else if(type == 4){
          
          Kmax<-length(f)
          T[i] <- t + sample(1:Kmax,1,prob = f)
          t<-T[i]
          i<- i+1
          
        }else{
          stop("The sojourn time is not valid")
        }
        
      }
      s<-i-1
      
      
    }else{
      ### param
      #fij
      if( type == 1 ){
        
        f<-param ## distribution parameters        
        if( !is.array(f) ){
          stop("f must be an array")
        }
        
        if( !is.matrix(distr) ){
          stop("The parameter distr must be a matrix")
        }
        
        while (t <= lengthSeq[m]){
          J[i + 1] <- sample(E, 1, prob = Ptrans[which(E==J[i]), ])
          #         cat("J")
          #         print(J[i])
          #         cat("T")
          #         print(t)
          
          
          if(distr[.code(J[i],E),.code(J[i+1],E)] == "unif"){
            ## Kmax = as.integer(1/f[.code(J[i],E), .code(J[i+1],E)])
            Kmax = f[.code(J[i],E),.code(J[i+1],E),1]
            T[i] <- t + sample(1:Kmax,1)
            t<-T[i]
            i<- i+1
            
            
          }else if(distr[.code(J[i],E),.code(J[i+1],E)] == "geom"){
            k<-rgeom(1,f[.code(J[i],E), .code(J[i+1],E),1])+1
            T[i] <- t + k
            t<-T[i]
            i<- i+1
            
          }else if(distr[.code(J[i],E),.code(J[i+1],E)] == "pois"){
            k<-rpois(1,f[.code(J[i],E), .code(J[i+1],E),1])+1 # no sojourn time equal to zero so we shift the Poisson distribution
            T[i] <- t + k
            t<-T[i]
            i<- i+1
            
          }else if(distr[.code(J[i],E),.code(J[i+1],E)] == "dweibull"){
            
            k<-rdweibull(1,f[.code(J[i],E), .code(J[i+1],E),1],f[.code(J[i],E), .code(J[i+1],E),2], zero = FALSE)
            #           cat("k weibull")
            #           print(k)
            T[i] <- t + k
            t<-T[i]
            i<- i+1
            
          }else if(distr[.code(J[i],E),.code(J[i+1],E)] == "nbinom"){
            k<-rnbinom(1,f[.code(J[i],E), .code(J[i+1],E),1], ,f[.code(J[i],E), .code(J[i+1],E),2])+1 ## we shift
            T[i] <- t + k
            t<-T[i]
            i<- i+1
            
          }else{
            stop("The distribution is not valid")
          }
          
        }
        #       s<-i-1
        #       JT<-rbind(J[1:s],T[1:s])
        #       y<-as.vector(.donneSeq(J,T))       
        #       return(y) 
        
        #fi
      }else if( type == 2 ){
        
        f<-param
        if( !is.matrix(f) ){
          stop("Param must be a matrix")
        }
        
        if( !is.vector(distr) ){
          stop("The distributions must be given in a vector")
        }
        
        while (t <= lengthSeq[m]){
          J[i + 1] <- sample(E, 1, prob = Ptrans[which(E==J[i]), ])
          
          if(distr[.code(J[i],E)] == "unif"){
            Kmax = f[.code(J[i],E),1] ##as.integer(1/f[.code(J[i],E)])
            T[i] <- t + sample(1:Kmax,1)#,prob = rep(f[.code(J[i],E)],Kmax))
            t <- T[i]
            i <- i+1
            
          }else if(distr[.code(J[i],E)] == "geom"){
            k<-rgeom(1,f[.code(J[i],E),1])+1 ## The 1 correspond to the first parameter and we shift
            T[i] <- t + k
            t<-T[i]
            i<- i+1
            
          }else if(distr[.code(J[i],E)] == "pois"){
            k<-rpois(1,f[.code(J[i],E),1])+1 ## The 1 correspond to the first parameter and we shift
            T[i] <- t + k
            t<-T[i]
            i<- i+1
            
          }else if(distr[.code(J[i],E)] == "dweibull"){
            k<-rdweibull(1,f[.code(J[i],E),1],f[.code(J[i],E),2], zero = FALSE)
            T[i] <- t + k
            t<-T[i]
            i<- i+1
            
          }else if(distr[.code(J[i],E)] == "nbinom"){
            k<-rnbinom(1,f[.code(J[i],E),1], ,f[.code(J[i],E),2])+1 ## The 1 correspond to the first parameter and we shift
            T[i] <- t + k
            t<-T[i]
            i<- i+1          
          }else{
            stop("The distribution is not valid")
          }
          
        }
        #       s<-i-1
        #       JT<-rbind(J[1:s],T[1:s])
        #       y<-as.vector(.donneSeq(J,T))       
        #       return(y) 
        
        #fj
      }else if( type == 3 ){
        f<-param
        if( !is.matrix(f) ){
          stop("param must be a matrix")
        }
        
        if( !is.vector(distr) ){
          stop("The distributions must be given in a vector")
        }
        
        while (t <= lengthSeq[m]){
          J[i + 1] <- sample(E, 1, prob = Ptrans[which(E==J[i]), ])
          
          if(distr[.code(J[i+1],E)] == "unif"){
            Kmax = f[.code(J[i+1],E),1] 
            T[i] <- t + sample(1:Kmax,1) 
            t <- T[i]
            i <- i+1
            
          }else if(distr[.code(J[i+1],E)] == "geom"){
            
            ## Change code ?
            k<-rgeom(1,f[.code(J[i+1],E),1])+1 ## we shift
            T[i] <- t + k
            t<-T[i]
            i<- i+1
            
          }else if(distr[.code(J[i+1],E)] == "pois"){
            
            ## Change code
            k<-rpois(1,f[.code(J[i+1],E),1])+1 ## we shift
            T[i] <- t + k
            t<-T[i]
            i<- i+1
            
          }else if(distr[.code(J[i+1],E)] == "dweibull"){
            
            k<-rdweibull(1,f[.code(J[i+1],E),1],f[.code(J[i+1],E),2], zero = FALSE)
            T[i] <- t + k
            t<-T[i]
            i<- i+1
            
          }else if(distr[.code(J[i+1],E)] == "nbinom"){
            k<-rnbinom(1,f[.code(J[i+1],E),1], ,f[.code(J[i+1],E),2])+1 ## we shift
            T[i] <- t + k
            t<-T[i]
            i<- i+1
            
          }else{
            stop("The distribution is not valid")
          }
          
        }
        #       s<-i-1
        #       JT<-rbind(J[1:s],T[1:s])
        #       y<-as.vector(.donneSeq(J,T))       
        #       return(y) 
        
        #f
      }else if( type == 4 ){
        
        f<-param
        if( is.list(f) || is.matrix(f) ){
          stop("param must be a vector or a number")
        }
        
        DISTR<-c("unif", "geom", "pois", "dweibull", "nbinom")
        distribution<-pmatch(distr, DISTR)
        
        while (t <= lengthSeq[m]){
          J[i + 1] <- sample(E, 1, prob = Ptrans[which(E==J[i]), ])
          
          ## unif
          if(distribution == 1){
            Kmax = f
            T[i] <- t + sample(1:Kmax,1)
            t <- T[i]
            i <- i+1
            
            ## geom
          }else if(distribution == 2){
            k<-rgeom(1,f)+1 ## we shift
            T[i] <- t + k
            t<-T[i]
            i<- i+1
            
            ## pois
          }else if(distribution == 3){
            k<-rpois(1,f)+1 ## we shift
            T[i] <- t + k
            t<-T[i]
            i<- i+1
            
            ## dweibull
          }else if(distribution == 4){
            k<-rdweibull(1,f[1],f[2], zero = FALSE)
            T[i] <- t + k
            t<-T[i]
            i<- i+1
            
            ## nbinom
          }else if(distribution == 5){
            k<-rnbinom(1,f[1], ,f[2])+1 ## we shift
            T[i] <- t + k
            t<-T[i]
            i<- i+1
            
          }else{
            stop("The distribution is not valid")
          }
          
        }
        s<-i-1
        
        #       cat("t")
        #       print(t)
        
        
        
      }else{
        stop("the Sojurn time is not valid")
      }
      
    }
    
    
    if(cens.beg == 1 && cens.end == 1){
      #     l = T[length(T)-1]-T[2]
      l = t -lengthSeq[m]
      #     cat("T")
      #     print(T)
      #     cat("l")
      #     print(l)
      n = lengthSeq[m]
      Nl = floor(l/2)
      #     cat("Nl")
      #     print(Nl)
      JT<-rbind(J[1:s],T[1:s])
      #     cat("seq")
      y<-as.vector(.donneSeq(J,T)) 
      #     print(length(y))
      y<-y[Nl:(t-1-Nl)]
      
      #First time is a Jump Time
    }else if(cens.beg == 0 && cens.end == 1){
      #   cat("Ok")
      #     cat("Ok")
      JT<-rbind(J[1:s],T[1:s])
      y<-as.vector(.donneSeq(J,T)) 
      y<-y[1:lengthSeq[m]]
      
    }else if(cens.beg == 1 && cens.end == 0){
      l = t -lengthSeq[m]
      JT<-rbind(J[1:s],T[1:s])
      y<-as.vector(.donneSeq(J,T)) 
      y<-y[l:(t-1)]
      
    }else{
      ## First and last times are jump times
      JT<-rbind(J[1:s],T[1:s])
      y<-as.vector(.donneSeq(J,T))   
    }
    
    if(!is.null(File.out)){
      if (file.exists(File.out)){
        write.fasta(sequences = y, names = "seq", file.out = File.out, open = "a")
      }else{
        write.fasta(y, names = "seq", file.out = File.out)
      }
    }
    
    out[[m]] = y
   
    }
  return(out) 
}
