
## __________________________________________________________
##
## .code
##
## __________________________________________________________
##

.code <- function(sequence,E)
{
  seq<-sequence
  s<-length(E)
  for (i in 1:s)
  {
    seq[seq==E[i]]<-i
  }
  seq<-as.numeric(seq)
  return(seq)
}


## __________________________________________________________
##
## .decode
##
## __________________________________________________________
##

.decode <-function(seq,E)
{
  sequence<-seq
  s<-length(seq)
  for (i in 1:s)
  {
    sequence[i]<-E[seq[i]]
  }
  return(sequence)    
}



## __________________________________________________________
##
## .donneJT
##
## __________________________________________________________
##

.donneJT <-function(seq)
{
  J<-NULL
  T<-NULL
  n<-length(seq)
  i1<-1
  i2<-2
  
  while (i1 <= n){
    while (i2<=n && (seq[i2]==seq[i1])){
      i2<-i2+1
    }
    J<-cbind(J,seq[i1])
    T<-cbind(T,i2)
    i1<-i2
    i2<-i1+1
  }
  JT <-rbind(J,T)
  colnames(JT)<-NULL
  return(JT)
}

## __________________________________________________________
##
## .donneSeq
##
## __________________________________________________________
##

.donneSeq <- function(J,T)
{
  seq<-NULL
  i<-1
  n<-length(T)
  t<-1
  while (i<= n)
  {
    k<-T[i]-t
    seq <-c(seq,rep(J[i],k))
    # for (j in 1:k)
    # {
    #   seq<-cbind(seq,J[i])
    # }
    t<-T[i]
    i<- i+1
  }
  
  return(seq)
}

## __________________________________________________________
##
## .donneYU
##
## __________________________________________________________
##

.donneYU<-function(JT, E){
  J<-.code(JT[1,],E)
  T<-as.numeric(JT[2,])
  L<-c(1,T)
  Y<-NULL
  U<-NULL
  n<-length(J)
  
  for (i in 1:n){
    for (j in 0:(L[i+1]-L[i]-1)){
      Y<-cbind(Y,J[i])
      U<-cbind(U,j)
    }
  }
  YU<-rbind(Y,U)
  colnames(YU)<-NULL
  return (YU)
}

## __________________________________________________________
##
## .comptage
##
## __________________________________________________________
##
.comptage <- function(J,L,S,Kmax){
  
  if(!is.list(J) || !is.list(L)){
    stop("J and L should be lists")
  }
  
  nbSeq = length(J)
  sNijk<-array(0,c(S,S,Kmax))
  Nstart = c()
  Nbij = array(0, c(S,S,Kmax))
  Nei = matrix(0, nrow = S, ncol = Kmax)
  Nbj = matrix(0, nrow = S, ncol = Kmax)
  Nb = vector( length = Kmax)
  Ne = vector( length = Kmax)  
  LsNijk<-list()
  LNstart = list()
  LNbij = list()
  LNei = list()
  LNbj = list()
  LNb = list()
  LNe = list()
  LNijk = list()
  for (k in 1:nbSeq){
    Ji = J[[k]]
    Li = L[[k]]
    M = length(Ji)
    Nstart = c(Nstart, Ji[1])
    Nbij[Ji[1],Ji[2],Li[1]] = Nbij[Ji[1],Ji[2],Li[1]] + 1
    LNbij[[k]] = Nbij
    Nei[Ji[M],Li[M]] = Nei[Ji[M],Li[M]] + 1
    LNei[[k]] = Nei
    # Nbj[Ji[2],Li[1]] = Nbj[Ji[2],Li[1]] + 1
    # Nb[Li[1]] = Nb[Li[1]] + 1
    # Ne[Li[M]] = Ne[Li[M]] + 1
    Nijk<-array(0,c(S,S,Kmax))
    for (i in 2:M){
      Nijk[Ji[i-1],Ji[i],Li[i-1]]<-Nijk[Ji[i-1],Ji[i],Li[i-1]]+1
    }
    sNijk <- sNijk + Nijk
    LNijk[[k]] = Nijk
  }
  
  LNij = list()
  Nij<-rowSums(sNijk,dims=2)
  for( k in 1:nbSeq ){
    LNij[[k]] = rowSums(LNijk[[k]], dims = 2)
  }
  
  diago = seq(1, S*S, by = S+1)
  Nij.temp = as.vector(Nij)[-diago]
  if( length(which(Nij.temp==0))!=0) {
    warning("Warning : missing transitions")
  }
  Nbj = colSums(Nbij)
  Ne = colSums(Nei)
  Nb = colSums(Nbj)
  Nbi = apply(Nbij, 3, rowSums)
  Nebi = Nei + Nbi
  Neb = Ne + Nb
  Nstarti = as.vector(count(seq = Nstart, wordsize = 1, alphabet = 1:S))
  LNstarti = list()
  for (k in 1:nbSeq){
    LNstarti[[k]] = as.vector(count(seq = Nstart[k], wordsize = 1, alphabet = 1:S))
  }
  Ni<-rowSums(Nij)
  if( length(which(Ni==0))!=0) {
    warning("Warning : missing letters")
  }
  Nj<-colSums(Nij)
  if( length(which(Nj==0))!=0) {
    warning("Warning : missing letters")
  }
  N<-sum(Nij)
  LNk = list()
  Nk<-colSums(sNijk,dims=2)
  for (k in 1:nbSeq){
    LNk[[k]] = colSums(LNijk[[k]], dims = 2)
  }
  LNik = list()
  Nik<-apply(sNijk,3,rowSums) #  i in raws and k in columns
  for (k in 1:nbSeq){
    LNik[[k]]<-apply(LNijk[[k]],3,rowSums)
  }
  LNjk = list()
  Njk<-colSums(sNijk) #  j in raws and k in columns
  for (k in 1:nbSeq){
    LNjk[[k]]<-colSums(LNijk[[k]])
  }
  return(list(Nijk=sNijk,Nij=Nij, Ni=Ni, Nj=Nj, N=N, Nk=Nk, Nik=Nik, Njk=Njk, Nstarti = Nstarti, Nbij = Nbij, Nei = Nei,
              Nbj = Nbj, Neb = Neb, Nebi = Nebi, Ne = Ne, Nb = Nb, Nbi = Nbi, LNijk = LNijk, LNij = LNij, LNij = LNij, LNik = LNik, LNjk = LNjk, 
              LNk = LNk, Nstart = LNstarti))
  
}

## __________________________________________________________
##
## .strationnary.law
##
## __________________________________________________________
##
# tpm: transition matrix p_ij = P(x_{t+1} = j | x_t = i)
.stationnary.law <- function(Pest){
  ## nb states
  m <- dim(Pest)[1]
  ## computation 
  A <- t(Pest) - diag(1, m, m)
  A[m, ] <- 1
  b <- c(rep(0, (m-1)), 1)
  stat.law <- solve(A, b)
  ## results
  return(stat.law)
}


## __________________________________________________________
##
## .limit.law
##
## __________________________________________________________
##
.limit.law <- function(q, pij){
  
  if ( dim(q)[1] != dim(pij)[1] && dim(pij)[2] != dim(q)[2] ){
    stop("The state number is not valid")
  }
  
  Kmax = dim(q)[3]
  S = dim(pij)[1]
  
  hj = apply(X = q, MARGIN = c(1,3), FUN = sum)
  ksup1 = rep(1:Kmax, S)
  mj_int = matrix(ksup1*as.vector(t(hj)), nrow = S, ncol = Kmax, byrow = TRUE)
  mj = rowSums(mj_int)
  statLaw = .stationnary.law(Pest = pij)
  
  Pi.j = statLaw*mj/sum(statLaw*mj)
  
  return(Pi.j)  
}


################################################################
## Computation of semi-Markov kernel (non-parametric case) ##
################################################################
## __________________________________________________________
##
## .calcul_q
##
## __________________________________________________________
##
.calcul_q = function(p,f,Kmax){
  
  dim.q=dim(f)
  q<-array(0,dim.q)
  for (k in 1:Kmax) {
    q[,,k] = p*f[,,k]
  }
  
  return(q)
  
}


############################################################
## Computation odf semi-Markovian kernel (parametric case)##
############################################################
### CASE fij ###
## __________________________________________________________
##
## .kernel_param_fij
##
## __________________________________________________________
##
# Computation of the kernel with the param just estimated
.kernel_param_fij <- function(distr, param, Kmax, pij, S){
  # get all the distributions in a vector
  lois = as.vector(distr)
  S2 = S*S
  fik = matrix(0, nrow = S2, ncol = Kmax)
  param = as.vector(param)
  diag = seq(from = 1, to = S*S, by = S+1)
  
  # print(S)
  # print(lois)
  # print(param)
  
  if ( "pois" %in% lois ){
    i = which(distr == "pois")
    for (j in i){
      fik[j,] = dpois(0:(Kmax-1), lambda = param[j])  
    }
  }
  if ( "geom" %in% lois ){
    i = which(distr == "geom")
    for (j in i){
      fik[j,] = dgeom(0:(Kmax-1), prob = param[j])
    }
  }
  if ( "nbinom" %in% lois ){
    i = which(distr == "nbinom")
    for (j in i){
      # print(j)
      # print(param[j])
      # print(param[j+S*S])
      fik[j,] = dnbinom(0:(Kmax-1), size = param[j], mu = param[j+S*S])
      if (fik[j,1] == 1){ fik[j,-1] = 0}else if (j %in% diag){ fik[j,] = 0}
      # print(fik)
      
    }
  }
  if ( "dweibull" %in% lois ){
    i = which(distr == "dweibull")
    for (j in i){
      fik[j,] = ddweibull(1:Kmax, q = param[j], beta = param[j+S*S])
      if (fik[j,1] == 1){ fik[j,-1] = 0}else if (j %in% diag){ fik[j,] = 0}
    }
  }
  fijk = array(fik, c(S,S,Kmax))
  q = array(pij, c(S,S,Kmax))*fijk
  
  return(list(q = q, f = fijk))
}

## __________________________________________________________
##
## .kernel_param_fj
##
## __________________________________________________________
##
### CASE fj ###
# Computation of the kernel with the param just estimated
.kernel_param_fj <- function(distr, param, Kmax, pij, S){
  
  lois = as.vector(distr)
  fik = matrix(data = 0, nrow = S, ncol = Kmax)
  
  if ( "pois" %in% lois ){
    i = which(distr == "pois")
    for (j in i){
      fik[j,] = dpois(0:(Kmax-1), lambda = param[j,1])
    }
    # i = which(distr == "pois")
    # fik[i,] = dpois(0:(Kmax-1), lambda = param[i])
  }
  if ( "geom" %in% lois ){
    i = which(distr == "geom")
    for (j in i){
      fik[j,] = dgeom(0:(Kmax-1), prob = param[j,1])
    }
  }
  if ( "nbinom" %in% lois ){
    i = which(distr == "nbinom")
    for (j in i){
      fik[j,] = dnbinom(0:(Kmax-1), size = param[j,1], mu = param[j,2])
    }
  }
  if ( "dweibull" %in% lois ){
    i = which(distr == "dweibull")
    for (j in i){
      fik[j,] = ddweibull(1:Kmax, q = param[j,1], beta = param[j,2])  
    }
  }
  
  f = rep(as.vector(t(fik)), each = S)
  fmat = matrix(f, nrow = Kmax, ncol =S*S, byrow = T)
  fk = array(as.vector(t(fmat)), c(S,S,Kmax))
  # fi.k = matrix(rep(as.vector(t(fik)), S), nrow = S*S, ncol = Kmax, byrow = TRUE)
  
  q = array(pij, c(S,S,Kmax))*fk
  q[which(is.na(q))]<-0
  
  return(list(q = q, f = fk))
}

## __________________________________________________________
##
## .kernel_param_fi
##
## __________________________________________________________
##
### CASE fi ###
# Computation of the kernel with the param just estimated
.kernel_param_fi <- function(distr, param, Kmax, pij, S){
  
  lois = as.vector(distr)
  fik = matrix(data = 0, nrow = S, ncol = Kmax)
  
  if ( "pois" %in% lois ){
    i = which(distr == "pois")
    for (j in i){
      fik[j,] = dpois(0:(Kmax-1), lambda = param[j,1])
    }
    # i = which(distr == "pois")
    # fik[i,] = dpois(0:(Kmax-1), lambda = param[i])
  }
  if ( "geom" %in% lois ){
    i = which(distr == "geom")
    for (j in i){
      fik[j,] = dgeom(0:(Kmax-1), prob = param[j,1])
    }
  }
  if ( "nbinom" %in% lois ){
    i = which(distr == "nbinom")
    for (j in i){
      fik[j,] = dnbinom(0:(Kmax-1), size = param[j,1], mu = param[j,2])
    }
  }
  if ( "dweibull" %in% lois ){
    i = which(distr == "dweibull")
    for (j in i){
      fik[j,] = ddweibull(1:Kmax, q = param[j,1], beta = param[j,2])  
    }
  }
  
  f = rep(as.vector(t(fik)), each = S)
  fmat = matrix(f, nrow = Kmax, ncol =S*S, byrow = TRUE)
  # fmat2 = matrix(as.vector(t(fmat)))
  fk = array(as.vector(t(fmat)), c(S,S,Kmax))
  fi.k = apply(X = fk,MARGIN =  c(1,3),FUN =  t)
  # fi.k = matrix(rep(as.vector(t(fik)), S), nrow = S*S, ncol = Kmax, byrow = TRUE)
  
  q = array(pij, c(S,S,Kmax))*fi.k
  q[which(is.na(q))]<-0
  
  return(list(q = q, f = fk))
}


## __________________________________________________________
##
## .kernel_param_f
##
## __________________________________________________________
##
### CASE f ###
.kernel_param_f <- function(distr, param, Kmax, pij, S){
  if ( distr == "pois" ) {
    f = dpois(0:(Kmax-1), lambda = param[1])
  }
  if ( distr == "geom" ) {
    f = dgeom(0:(Kmax-1), prob = param[1])
  }
  if ( distr == "dweibull" ) {
    f = ddweibull(1:Kmax, q = param[1], beta = param[2])
  }
  if ( distr == "nbinom" ){
    f = dnbinom(0:(Kmax-1), size = param[1], mu = param[2])
  }
  
  f = rep(f, each = S*S)
  fmat = matrix(f, nrow = Kmax, ncol =S*S, byrow = TRUE)
  fk = array(as.vector(t(fmat)), c(S,S,Kmax))
  # fk = array(fmat, c(S,S,Kmax)) 
  q = array(pij, c(S,S,Kmax))*fk
  q[which(is.na(q))]<-0
  
  return(list(q = q, f = fk))
}



##################################
## NON-PARAMETRIC (couple MC) ##
##################################
## __________________________________________________________
##
## .calcul.Niujv.matrice 
##
## __________________________________________________________
##
##
.calcul.Niujv.matrice <- function(Y,U,S,Kmax){
  
  if(!is.list(Y) || !is.list(U)){
    stop("Y and U should be lists")
  }
  
  nbSeq = length(Y)
  sNiujv<-matrix(0, nrow=S*Kmax, ncol=S*Kmax)
  for (k in 1:nbSeq){
    Yi = Y[[k]]
    Ui = U[[k]]
    M = length(Yi)
    Niujv=matrix(0, nrow=S*Kmax, ncol=S*Kmax)
    for (i in 2:M){
      if((Yi[i-1] != Yi[i] && Ui[i] == 0) || (Yi[i-1] == Yi[i] && Ui[i] == Ui[i-1] + 1)){
        Niujv[Yi[i-1]+(Yi[i-1]-1)*(Kmax-1)+Ui[i-1],Yi[i]+(Yi[i]-1)*(Kmax-1)+Ui[i]]<-
          Niujv[Yi[i-1]+(Yi[i-1]-1)*(Kmax-1)+Ui[i-1],Yi[i]+(Yi[i]-1)*(Kmax-1)+Ui[i]]+1
      }
    }
    
    #     print(Niujv)
    sNiujv <- sNiujv + Niujv
  }
  
  return(sNiujv)
  
}

## __________________________________________________________
##
## .calcul.Niujv 
##
## __________________________________________________________
##
.calcul.Niujv <- function(Y,U,S,Kmax){
  
  if(!is.list(Y) || !is.list(U)){
    stop("Y and U should be lists")
  }
  
  nbSeq = length(Y)
  for (k in 1:nbSeq){
    ## Shift from aphabet {1,2,3,4,...} to {0,1,2,3,...} 
    Yi = Y[[k]]-1
    Ui = U[[k]]
    M = length(Yi)
    Niujv=rep(0, (Kmax+1)*S*(S-1)+S*(Kmax))
    for (i in 2:M){
#       if((Yi[i-1] != Yi[i] && Ui[i] == 0) || (Yi[i-1] == Yi[i] && Ui[i] == Ui[i-1] + 1)){
#         Niujv[Yi[i-1]+(Yi[i-1]-1)*(Kmax-1)+Ui[i-1],Yi[i]+(Yi[i]-1)*(Kmax-1)+Ui[i]]<-
#           Niujv[Yi[i-1]+(Yi[i-1]-1)*(Kmax-1)+Ui[i-1],Yi[i]+(Yi[i]-1)*(Kmax-1)+Ui[i]]+1
#       }
      if( Yi[i-1] != Yi[i] ){
        if( Ui[i] != 0 ){
          stop("erreur")
        }
        if( Yi[i] > Yi[i-1] ){
          Niujv[Yi[i-1]*(Kmax+1)*(S-1)+(Kmax+1)*(Yi[i]-1)+Ui[i-1]+1] =
            Niujv[Yi[i-1]*(Kmax+1)*(S-1)+(Kmax+1)*(Yi[i]-1)+Ui[i-1]+1] + 1
        }else{
          Niujv[Yi[i-1]*(Kmax+1)*(S-1)+(Kmax+1)*Yi[i]+Ui[i-1]+1] = 
            Niujv[Yi[i-1]*(Kmax+1)*(S-1)+(Kmax+1)*Yi[i]+Ui[i-1]+1] + 1
        } 
      }else{
        if( Ui[i] != Ui[i-1] + 1 ){
          stop("erreur")
        }
        Niujv[(Kmax+1)*S*(S-1)+Yi[i-1]*Kmax+Ui[i-1]+1] = Niujv[(Kmax+1)*S*(S-1)+Yi[i-1]*Kmax+Ui[i-1]+1] + 1 
      }
    }
    
    #     print(Niujv)
#     sNiujv <- sNiujv + Niujv
  }
  
  return(Niujv)
  
}

## __________________________________________________________
##
## .calcul.Niu
##
## __________________________________________________________
##
.calcul.Niu = function(Kmax,S,Niujv){

  Niu<-matrix(0, nrow=S, ncol=Kmax+1)
  for (i in 0:(S-1)){
    for (u in 0:(Kmax-1)){
      Niu[i+1,u+1] = sum(Niujv[i*(Kmax+1)*(S-1)+(Kmax+1)*(0:(S-2))+u+1]) + 
        Niujv[(Kmax+1)*S*(S-1)+i*Kmax+u+1]
          
    }
    Niu[i+1,Kmax] = sum(Niujv[i*(Kmax+1)*(S-1)+(Kmax+1)*(0:(S-2))+u+1])
  }  
    
  return (Niu)
  
}

## __________________________________________________________
##
## .p.hat 
##
## __________________________________________________________
##
.p.hat = function(Kmax,S,Niujv,Niu){
  
  p = rep(0, (Kmax+1)*S*(S-1)+S*(Kmax))
  for (i in 0:(S-1)){
    for (u in 0:(Kmax-1)){
      for (j in 0:(S-2)){
        #         for (v in 0:(Kmax-1)){
        if( i != j ){
          if( i < j ){
            if( Niu[i+1,u+1] == 0 ){
              break
            }
            p[(Kmax+1)*i*(S-1)+(j-1)*(Kmax+1)+u+1] = 
              Niujv[(Kmax+1)*i*(S-1)+(j-1)*(Kmax+1)+u+1]/Niu[i+1,u+1]
            
          }else{
            if( Niu[i+1,u+1] == 0 ){
              break
            }
            p[(Kmax+1)*i*(S-1)+j*(Kmax+1)+u+1] = 
              Niujv[(Kmax+1)*i*(S-1)+j*(Kmax+1)+u+1]/Niu[i+1,u+1]
            
          }
        }else{
          if( Niu[i+1,u+1] == 0 ){
            break
          }
          #             if( v == u + 1 ){
          
          p[(Kmax+1)*S*(S-1)+i*Kmax+u+1] = Niujv[(Kmax+1)*S*(S-1)+i*Kmax+u+1]/Niu[i+1,u+1]
          
          #             }
        }
        
        #           }
      }
    }      
    
  }
  return (p)
}  


## __________________________________________________________
##
## .calcul_NP.q
##
## __________________________________________________________
##
.calcul_NP.q = function(S,Kmax,p){
  
  q = array(0, c(S,S,Kmax))
  Pi = 1
  for (i in 0:(S-1)){
    for (j in 0:(S-1)){
      if (i!=j){
      for (k in 1:Kmax){
        Pr = 1        
       
        if ( i > j ){
          Pi = p[(Kmax+1)*i*(S-1)+j*(Kmax+1)+k] #(k-1)+1
          #           cat("Pi :")
          #           print(Pi)
          #           cat("Pr :")
          #           print(Pr)
          #           cat("==================================")
          #           cat("\n")

        }
        if ( i < j ){
          Pi = p[(Kmax+1)*i*(S-1)+(j-1)*(Kmax+1)+k] #(k-1)+1
          #           cat("Pi :")
          #           print(Pi)
          #           cat("Pr :")
          #           print(Pr)
          #           cat("==================================")
          #           cat("\n")

        }
        Pr=Pi*Pr
        if( k >= 2 ){
          for (t in 0:(k-2)){
            
            Pr = Pr * p[(Kmax+1)*S*(S-1)+ i*Kmax + t+1]
            
          }
        }
        
        q[i+1,j+1,k]=Pr
        
      }
      }
    }
  }
  return (q)
}

## __________________________________________________________
##
## .calcul.NP.q.matrice
##
## __________________________________________________________
##
.calcul_NP.q.matrice = function(Y,U,p,S,Kmax){
  if(!is.list(Y) || !is.list(U)){
    stop("Y and U should be lists")
  }
  
  nbSeq = length(Y) 
  q=array(0, c(S,S,Kmax))
  for (l in 1:nbSeq){
    Yi = Y[[l]]
    Ui = U[[l]]
    M = length(Yi)
    
    for (i in 2:M){
      if((Yi[i-1] != Yi[i] && Ui[i] == 0) || (Yi[i-1] == Yi[i] && Ui[i] == Ui[i-1] + 1)){
        for (k in 1:Kmax){
          Pr=1
          if( k >= 2){
            for (j in 1:(k-1)){
              Pr = Pr*p[Yi[i-1]+(Yi[i-1]-1)*(Kmax-1)+j-1,Yi[i-1]+(Yi[i-1]-1)*(Kmax-1)+j]
            }
          }
          q[Yi[i-1],Yi[i],k]=Pr*p[Yi[i-1]+(Yi[i-1]-1)*(Kmax-1)+k-1,Yi[i]+(Yi[i]-1)*(Kmax-1)]
        }
      }
    }
    
  }
  return (q)
}

## __________________________________________________________
##
## .write.results
##
## __________________________________________________________
##
.write.results = function(TypeSojournTime, cens.beg, cens.end, p, param){
  date.res = paste("results", Sys.Date(), sep="_")
  File.out = paste(date.res, "txt", sep=".")
  i = 1
  while(file.exists(File.out)){
    # print(File.out)
    ## get the file name before extension
    File.txt = unlist(strsplit(File.out, split = "[.]"))[1]
    if (i>1) {
      ## get the file name without version number
      File.1 = paste(unlist(strsplit(File.txt, split = "_"))[1], unlist(strsplit(File.txt, split = "_"))[2], sep="_")
      Filei = paste(File.1, i, sep="_")
      # Filei = paste(File.int, "txt", sep=".")
    }else{
      Filei = paste(File.txt, i, sep="_")
    }
    File.out = paste(Filei, "txt", sep=".")
    i = i + 1
  }
  sink(File.out)
  cat("The Results of estimation for one or several sequences \n")
  cat("\n")
  cat("The sojourn time chosen : ")
  cat("\n")
  print(TypeSojournTime)
  # cat("The matrix of the distributions")
  # cat("\n")
  # print(distr)
  cat("Censoring at the beginnig")
  cat("\n")
  if(cens.beg == 1){print("yes")}else{print("no")}
  cat("Censoring at the end")
  cat("\n")
  if(cens.end == 1){print("yes")}else{print("no")}
  cat("The transition probabilities : ")
  cat("\n")
  print(p)
  cat("The estimated parameters : ")
  cat("\n")
  print(param)
  sink()
} 


## __________________________________________________________
##
## .estimNP_aucuneCensure
##
## __________________________________________________________
##
####################
##### INPUT: 
## seq: list of sequences
## E: state space
## TypeSojournTime: Type of Sojourn time ("fij","fi","fj","f")
## distr = "NP"
## No censoring
## All trajectories with the same type of censoring
#####################
##### OUTPUT: 
## p: transition matrix
## f: array
## q: array
.estimNP_aucuneCensure<-function(file = NULL, seq, E, TypeSojournTime="fij", distr="NP"){
  
  
  if(is.null(file) && is.null(seq)){
    stop("One of the two parameters file or seq should be given")
  }
  
  if(!is.null(file)){
    fasta = read.fasta(file = file)
    seq = getSequence(fasta)
    attr = getAnnot(fasta)
  }
  
  
  if(!is.list(seq)){
    stop("The parameter seq should be a list")
  }
  
  ## States
  TYPES<-c("fij", "fi", "fj", "f")
  type<-pmatch(TypeSojournTime,TYPES)
  
  ## Length of the state space
  S<-length(E)
  
  ## All sequences in a vector
  nbSeq<-length(seq)
  #  seq1<-c()
  #  for (k in 1:nbSeq){
  #    seq1<-c(seq1, unlist(seq[[k]]))
  #  }
  
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
    L[[k]]<-T[[k]] - c(1,T[[k]][-length(T[[k]]-1)]) ## sojourn time
    Kmax <- max(Kmax, max(L[[k]])) 
  }
  
  
  ## J: list
  ## L: list
  ## S: state space size
  ## Kmax: maximal sojourn time
  ################
  ## get the counts
  res<-.comptage(J,L,S,Kmax)
  
  Nij = res$Nij
  Ni = res$Ni
  Nj = res$Nj
  N= res$N
  Nk = res$Nk
  Nik = res$Nik
  Njk = res$Njk
  Nijk = res$Nijk
  Nstart = res$Nstarti
  Nbij = res$Nbij
  Nei = res$Nei
  Nbj = res$Nbj
  Neb = res$Neb
  Nebi = res$Nebi
  Nbi = res$Nbi
  Nb = res$Nb
  Ne = res$Ne
  #########
  
  
  
  if(is.na(distr)){
    stop("invalid distribution")
    
    ## Non-parametric
  }else if(distr == "NP"){
    
    if( type == 4 ){
      # Nij<-rowSums(Nijk,dims=2)
      # Ni<-rowSums(Nij)
      # N<-sum(Nij)
      mat.Nk<-matrix(1,nrow=S,ncol=S)
      diag(mat.Nk)<-0
      Nk<-array(mat.Nk,c(S,S,Kmax))
      
      ## Estimators
      if(length(which(Ni==0))!=0){ # Presence of all states
        
        print("Warning : All the states are not observed")
        cat("\n")
        
        p<-Nij/(as.vector(Ni))
        p[which(is.na(p))]<-0
        
        f <- Nk/N
        f[which(is.na(f))]<-0
        
        q<-Nijk/array(Ni,c(S,S,Kmax))
        q[which(is.na(q))]<-0
        
      } else {
        diago = seq(1, S*S, by = S+1)
        Nij.temp = as.vector(Nij)[-diago]
        if( length(which(Nij.temp==0))!=0) # Presence of all transitions
        {
          print("Warning: The transitions are not all observed")
          cat("\n")
          
          p<-Nij/(as.vector(Ni))
          
          f <- Nk/N
          f[which(is.na(f))]<-0
          
          q<-Nijk/array(Ni,c(S,S,Kmax))                    
        } else {
          p<-Nij/(as.vector(Ni))
          f <- Nk/N
          q<-Nijk/array(Ni,c(S,S,Kmax))
        }
      }
      # return (list(p,f,q))
      
    }else {
      
      if( type == 2 ){
        # Nij<-rowSums(Nijk,dims=2)
        Nik<-vector("list", Kmax)
        N<-vector("list", Kmax)
        list.Nijk<-lapply(seq(dim(Nijk)[3]), function(x) Nijk[ , , x])
        for (k in 1:Kmax){
          N[[k]]<-apply(list.Nijk[[k]],1,sum)
          Nik[[k]]<-matrix(unlist(N[[k]]),4,4)
          diag(Nik[[k]])<-0
        }
        Nik<-array(unlist(Nik), c(4,4,Kmax))
        Ni<-rowSums(Nij)
        
        ## Estimators
        if(length(which(Ni==0))!=0){ # Presence of all states
          
          print("Warning: Not all the states are observed")
          cat("\n")
          
          p<-Nij/(as.vector(Ni))
          p[which(is.na(p))]<-0
          
          f<-Nik/Ni
          f[which(is.na(f))]<-0
          
          q<-Nijk/array(Ni,c(S,S,Kmax))
          q[which(is.na(q))]<-0
          
        } else {
          diago = seq(1, S*S, by = S+1)
          Nij.temp = as.vector(Nij)[-diago]
          if( length(which(Nij.temp==0))!=0) # Presence of all transitions
          {
            print("Warning: The transitions are not all observed")
            cat("\n")
            
            p<-Nij/(as.vector(Ni))
            
            f<-Nik/Ni
            f[which(is.na(f))]<-0
            
            q<-Nijk/array(Ni,c(S,S,Kmax))                    
          } else {
            p<-Nij/(as.vector(Ni))
            f<-Nik/Ni
            q<-Nijk/array(Ni,c(S,S,Kmax))
          }
        }
        # return (list(p,f,q))
      }
      
      if( type == 1 ){
        # Nij<-rowSums(Nijk,dims=2)
        # Ni<-rowSums(Nij)
        ## Estimators
        
        if(length(which(Ni==0))!=0){ # Presence of all states
          
          print("Warning: The transitions are not all observed")
          cat("\n")
          
          p<-Nij/(as.vector(Ni))
          p[which(is.na(p))]<-0
          
          f <- Nijk/array(Nij,c(S,S,Kmax))
          f[which(is.na(f))]<-0
          
          q<-Nijk/array(Ni,c(S,S,Kmax))
          q[which(is.na(q))]<-0
          
        } else {
          diago = seq(1, S*S, by = S+1)
          Nij.temp = as.vector(Nij)[-diago]
          if( length(which(Nij.temp==0))!=0) # Presence of all transitions
          {
            print("Warning: The transitions are not all observed")
            cat("\n")
            
            p<-Nij/(as.vector(Ni))
            
            f <- Nijk/array(Nij,c(S,S,Kmax))
            f[which(is.na(f))]<-0
            
            q<-Nijk/array(Ni,c(S,S,Kmax))                    
          } else {
            p<-Nij/(as.vector(Ni))
            f <- Nijk/array(Nij,c(S,S,Kmax))
            q<-Nijk/array(Ni,c(S,S,Kmax))
          }
        }
        # return(list(p,f,q))
      }
      
      if( type == 3 ){   
        # Nij<-rowSums(Nijk,dims=2)
        # Ni<-rowSums(Nij)
        Njk<-vector("list", Kmax)
        N<-vector("list", Kmax)
        list.Nijk<-lapply(seq(dim(Nijk)[3]), function(x) Nijk[ , , x])
        for (k in 1:Kmax){
          N[[k]]<-apply(list.Nijk[[k]],2,sum)
          Njk[[k]]<-matrix(unlist(N[[k]]),4,4)
          diag(Njk[[k]])<-0
        }
        Njk<-array(unlist(Njk), c(4,4,Kmax))
        Nj<-colSums(Nij)
        ## Estimators
        
        
        if(length(which(Ni==0))!=0){ # Presence of all states
          
          print("Warning: The transitions are not all observed")
          cat("\n")
          
          p<-Nij/(as.vector(Ni))
          p[which(is.na(p))]<-0
          
          f<-Njk/Nj
          f<-lapply(seq(dim(f)[3]), function(x) f[ , , x])
          Lf<-vector("list",5)
          for (k in 1:Kmax){
            Lf[[k]]<-matrix(f[[k]],4,4)
            f[[k]]<-t(Lf[[k]])
          }
          f<-array(unlist(f), c(4,4,Kmax))
          f[which(is.na(f))]<-0
          
          q<-Nijk/array(Ni,c(S,S,Kmax))
          q[which(is.na(q))]<-0
          
        } else {
          diago = seq(1, S*S, by = S+1)
          Nij.temp = as.vector(Nij)[-diago]
          if( length(which(Nij.temp==0))!=0) # Presence of all transitions
          {
            print("Warning: The transitions are not all observed")
            cat("\n")
            
            p<-Nij/(as.vector(Ni))
            
            f<-Njk/Nj
            f<-lapply(seq(dim(f)[3]), function(x) f[ , , x])
            Lf<-vector("list",5)
            for (k in 1:Kmax){
              Lf[[k]]<-matrix(f[[k]],4,4)
              f[[k]]<-t(Lf[[k]])
            }
            f<-array(unlist(f), c(4,4,Kmax))
            f[which(is.na(f))]<-0
            
            q<-Nijk/array(Ni,c(S,S,Kmax))                    
          } else {
            p<-Nij/(as.vector(Ni))
            f<-Njk/Nj
            f<-lapply(seq(dim(f)[3]), function(x) f[ , , x])
            Lf<-vector("list",5)
            for (k in 1:Kmax){
              Lf[[k]]<-matrix(f[[k]],4,4)
              f[[k]]<-t(Lf[[k]])
            }
            f<-array(unlist(f), c(4,4,Kmax))
            q<-Nijk/array(Ni,c(S,S,Kmax))
          }
          
        }
        
    
        
      }
      
    }
    
  }
  
  ## Inital law
  if(nbSeq > 1){
    init1 = Nstart/sum(Nstart)
  }else{
    ## Computation of the limit distribution
    init1 = .limit.law(q = q, pij = p)
  }
  
  f[which(is.na(f))]<-0
  return(list(init = init1, Ptrans = p, laws = f, q = q))
  
}

## __________________________________________________________
##
## .estimNP_CensureFin
##
## __________________________________________________________
##
## Estimation with the associated MC (Y,U) for several trajectories
.estimNP_CensureFin<-function(file = NULL, seq, E, TypeSojournTime){

  
  if(is.null(file) && is.null(seq)){
    stop("One of the two parameters file or seq must be given")
  }
  
  if(!is.null(file)){
    fasta = read.fasta(file = file)
    seq = getSequence(fasta)
    attr = getAnnot(fasta)
  }
  
  
  if(!is.list(seq)){
    stop("The parameter seq must be a list")
  }
  
  ## States
  TYPES<-c("fij", "fi", "fj", "f")
  type<-pmatch(TypeSojournTime,TYPES)
  
  ## Size of the state space
  S<-length(E)
  
  ## All sequences 
  nbSeq<-length(seq)
  KmaxStart = 0
  Kmax = 0
  Y<-list()
  U<-list()
  for (k in 1:nbSeq){
    ## Manipulations
    JT<-.donneJT(unlist(seq[[k]]))
    J[[k]]<-.code(JT[1,],E)
    T[[k]]<-as.numeric(JT[2,])
    L[[k]]<-T[[k]] - c(1,T[[k]][-length(T[[k]]-1)]) ## sojourn times
    KmaxStart <- max(KmaxStart, max(L[[k]])) 
    YU<-.donneYU(JT, E)
    Y[[k]]<-YU[1,]
    U[[k]]<-YU[2,]
    Kmax <- max(Kmax, max(U[[k]])+1) 
  }
  
  ## get the counts in first position for each state 
  res<-.comptage(J,L,S,KmaxStart)
  Nstart = res$Nstarti
  
  ## Computation of Niujv
  ## with array ?
  Niujv = .calcul.Niujv(Y,U,S,Kmax)
  
  #   Niujv.mat = .calcul.Niujv.matrice(Y,U,S,Kmax)
  
  #   Niu.mat<-rowSums(Niujv.mat)
  Niu = .calcul.Niu(Kmax,S,Niujv)
  
  #   p.mat<-Niujv.mat/(as.vector(Niu.mat))
  #   p.mat[which(is.na(p.mat))]<-0
  p = .p.hat(Kmax = Kmax, S = S, Niujv = Niujv, Niu = Niu)
  
  ##computation of qs  
  ## q.old = .calcul_NP.q.matrice(Y,U,p.mat,S,Kmax)
  q = .calcul_NP.q(S = S, Kmax = Kmax, p = p)
  
  
  pij = rowSums(q,dims=2)

  # fij
  if( type == 1 ){
    f=q/array(pij, c(S,S,Kmax))
    f[which(is.na(f))]<-0   
    
    # fi
  }else if( type == 2 ){
    # sum over the raws of each matrix of the array
    a = apply(q, c(1,3), sum)
    # transformation of the matrix to an array
    f = array(apply(a,2,function(x) matrix(x, nrow = S, ncol = S)), c(S,S,Kmax))
    # fj
  }else if( type == 3 ){
    a = apply(q, c(2,3),sum)
    sum_qi = array(0, c(S,S,Kmax))
    sum_qi = array(apply(a,2,function(x) matrix(x, nrow = S, ncol = S)), c(S,S,Kmax))
    mat_p = colSums(pij)
    calcul.f = array(apply(sum_qi,c(2,3),function(x) x/mat_p), c(S,S,Kmax))
    calcul.f[which(is.na(calcul.f))]<-0
    f = array(apply(calcul.f, 3, function(x) t(x)), c(S,S,Kmax))
    #     for (i in 1:Kmax){
    #       for( j in 1:S ){
    #         sum_qi[,,i] = matrix(colSums(q)[,i])
    #         sum_qi[,,i] = t(sum_qi[,,i])
    #         f[,j,i]  = sum_qi[,j,i]/colSums(p)[j]
    #       }
    #     }
    
    # f
  }else if ( type == 4 ){
    ## warning S = 3 (The states are not all observed)
    sum_qij = apply(q,3,sum)/S
    f = array(apply(as.matrix(sum_qij),1,function(x) matrix(x, nrow = S, ncol = S)),
              c(S,S,Kmax))
    
  }else{
    stop("The parameter \"TypeSojournTime\" is not valid")
  }
  
  ##Initial law
  if(nbSeq > 1){
    init = Nstart/sum(Nstart)
  }else{
    ## Computation of limit law
    init = .limit.law(q = q, pij = pij)
  }
  
  return (list(init = init, Ptrans = pij, laws = f, q = q))
  
}

## __________________________________________________________
##
## .estim.plusTraj
##
## __________________________________________________________
##
####################
##### INPUT: 
## file: path to the file
## seq: list of sequences
## E: state space
## TypeSojournTime: Type of Sojour time ("fij","fi","fj","f")
## distr:  - matrix of distributions if fij
##         - vector of distributions if  fi or fj
## cens.end = 1 if censoring at the end
##            0 if not
## cens.beg = 1 if censoring at the beginning 
##            0 if not
#####################
##### OUTPUT : 
## p: transition matrix
## param:  - array (in param[,,1] the first parameter of the distribution and in param[,,2] the second) if fij
##         - matrix (in param[,1] the first parameter of the distribution and in param[,2] the second) if fi ou fj
##         - vector (in param[1] the first parameter of the distribution and in param[2] the second) if f
.estim.plusTraj<-function(seq, E, TypeSojournTime = "fij", distr = "NP", cens.end = 0, cens.beg = 0){
  
  # require(seqinr)
  # require(DiscreteWeibull)
  # 

  
  if(!is.list(seq)){
    stop("the parameter seq should be a list")
  }
  # 
  ## States
  TYPES<-c("fij", "fi", "fj", "f")
  type<-pmatch(TypeSojournTime,TYPES)
  
  ## Length of state space
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
  ## S: alphabet size
  ## Kmax: maximal sojourn time
  ################
  ## get the counts
  res<-.comptage(J,L,S,Kmax)
  
  Nij = res$Nij
  Ni = res$Ni
  Nj = res$Nj
  N= res$N
  Nk = res$Nk
  Nik = res$Nik
  Njk = res$Njk
  Nijk = res$Nijk
  Nstart = res$Nstarti
  Nbij = res$Nbij
  Nei = res$Nei
  Nbj = res$Nbj
  Neb = res$Neb
  Nebi = res$Nebi
  Nbi = res$Nbi
  Nb = res$Nb
  Ne = res$Ne
  #########
  
  
  ## Initialization of P
  P <- matrix(0,S,S)
  
  ## fij
  if ( type == 1 ){
    
    if ( !is.matrix(distr) ){
      stop("The parameter distr should be a matrix")
    }
    
    param = array(0, c(S,S,2))
    ## Warning unif infinite ?
    for (i in 1:S){
      for ( j in 1:S ){
        
        if(i != j){
           if ( distr[i,j] == "unif" ){
             
            ## No censoring
            if( cens.end == 0 && cens.beg == 0 ){
              P = Nij/Ni
              
              maxUnif = which.max(Nijk[i,j,])
              param[i,j,1] = maxUnif
              
              ## Censoring at the end
            }else if ( cens.end == 1 && cens.beg == 0 ){
              ## Particular case S = 2 because index 0 for P
              if ( S == 2 ){
                #############################
                ## Estimation
                #############################
                plus = seq(from = 0, to = Kmax*S*S-1, by = S*S)
                ## vect.Nijk is the vector of the array Nijk
                vect.Nijk = as.vector(Nijk)
                ## the Nijk[i,j,] correspond for the vector Nijk to vect.Nijk[i+S*(j-1)+plus]
                ## estimation of p_{jt+1v} et theta_{jt+1v}
                Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
                El.diag = seq(from= 1, to = S*S, by = S+1)
                Nij.vect = Nij.vect0[-El.diag]
                # Nij0 = Nij0[c(etat.cens),]
                
                
                ### Separation of parameters in order to compute the log-likeihood
                Nij1 = Nij.vect[1]
                Nij2 = Nij.vect[2]
                lois = distr[i,]
                lois = lois[-i]
                
                ## Function to estimate
                logvrais = function(par){
                  p = NULL
                  g = NULL
                  nb = NULL
                  dw = NULL
                  u = NULL
                  
                  if("unif" %in% lois){
                    u = 2*which(lois == "unif")-1
                    u = c(u, u+1)
                  }
                  if("pois" %in% lois){
                    p = 2*which(lois == "pois")-1
                    p = c(p, p+1)
                  }
                  if("geom" %in% lois){
                    g = 2*which(lois == "geom")-1
                    g = c(g, g+1)
                  }
                  if("nbinom" %in% lois){
                    nb = 2*which(lois == "nbinom")-1
                    nb = c(nb, nb+1)
                  }
                  if("dweibull" %in% lois){
                    dw = 2*which(lois == "dweibull")-1
                    dw = c(dw, dw+1)
                  }
                  
                  if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                  
                  MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])]), 
                                 dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])]),
                                 dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])]), 
                                 ddw,
                                 dunif(1:Kmax, min = 0, max = par[u[1]])), nrow = Kmax)
                  
                  
                  # print(MGM)
                  fuv = matrix(nrow = S-1, ncol = Kmax)
                  if("pois" %in% lois){
                    nr = which(lois == "pois")
                    fuv[nr,] = MGM[,1]
                  }
                  if("geom" %in% lois){
                    nr = which(lois == "geom")
                    fuv[nr,] = MGM[,2]
                  }
                  if("nbinom" %in% lois){
                    nr = which(lois == "nbinom")
                    fuv[nr,] = MGM[,3]
                  }
                  if("dweibull" %in% lois){    
                    nr = which(lois == "dweibull")
                    fuv[nr,] = MGM[,4]
                  } 
                  if("unif" %in% lois){    
                    nr = which(lois == "unif")
                    fuv[nr,] = MGM[,5]
                  }
                  
                  puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                  puv.vect0 = as.vector(t(puv))
                  puv.vect = puv.vect0[-El.diag]
                  
                  -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                     # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                     +  sum(vect.Nijk[i+S*(j-1)+plus]*dunif(1:Kmax, min = 0, max = par[u[1]], log = TRUE))
                     + sum( Nei[i,]* log(1-puv.vect%*%fuv))
                  )  
                  
                  
                }
                
                ### Vector for the initialisation of optim
                init = c()
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                  init[p] = c(2, 1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                  init[g] = c(0.5, 0.4)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                  init[nb] = c(2, 4)
                }
                if("dweibull" %in% lois){ 
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                  init[dw] = c(0.6, 0.8)
                }
                if("unif" %in% lois){ 
                  u = 2*which(lois == "unif")-1
                  u = c(u, u+1)
                  init[u] = c(Kmax, Kmax)
                }
                
                
                ### Likelihood Optimisation
                CO2<-optim(par = init,fn = logvrais, method ="Nelder-Mead")
                ### Get parameters
                parC = CO2$par
                maxUnif = parC[u[1]]
                param[i,j,1] = maxUnif
                
                puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                P = puv
                
              }else{
                #############################
                ## Estimation
                #############################
                plus = seq(from = 0, to = Kmax*S*S-1, by = S*S)
                ## vect.Nijk est le vecteur de l'array Nijk
                vect.Nijk = as.vector(Nijk)
                ##  Nijk[i,j,] corresponds to the Nijk vector for vect.Nijk[i+S*(j-1)+plus]
                ## estimation function of  p_{jt+1v} et theta_{jt+1v}
                Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
                El.diag = seq(from= 1, to = S*S, by = S+1)
                Nij.vect = Nij.vect0[-El.diag]
                # Nij0 = Nij0[c(etat.cens),]
                
                
                ### Separation of parameters in order to compute log-likeihood
                Nij1 = c()
                Nij2 = c()
                a = 1
                while (a < length(Nij.vect)){
                  Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## corresponds to the (S-2) first values of Nij on each line
                  Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## corresponds to the last value of Nij on each line
                  a = a + (S-1)
                }
                
                lois = distr[i,]
                lois = lois[-i]
                
                ## Function to estimate
                logvrais = function(par){
                  m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                  p = cbind(m, 1-rowSums(m))
                  puv = p[i,]
                  # ind1 = seq(from = (S*(S-2))+1, to =(S*(S-2)+5), by=2)
                  # ind2 = seq(from = (S*(S-2))+2, to =(S*(S-2)+6), by=2)
                  p = NULL
                  g = NULL
                  nb = NULL
                  dw = NULL
                  u = NULL
                  
                  if("pois" %in% lois){
                    p = 2*which(lois == "pois")-1
                    p = c(p, p+1)
                  }
                  if("geom" %in% lois){
                    g = 2*which(lois == "geom")-1
                    g = c(g, g+1)
                  }
                  if("nbinom" %in% lois){
                    nb = 2*which(lois == "nbinom")-1
                    nb = c(nb, nb+1)
                  }
                  if("dweibull" %in% lois){
                    dw = 2*which(lois == "dweibull")-1
                    dw = c(dw, dw+1)
                  }
                  if("unif" %in% lois){
                    u = 2*which(lois == "unif")-1
                    u = c(u, u+1)
                  }
                  
                  if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                  
                  MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])]), 
                                 dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])]),
                                 dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])]), 
                                 ddw, 
                                 dunif(1:Kmax, min = 0, max = par[S*(S-2)+u[1]])), nrow = Kmax)
                  
                  # print(MGM)
                  fuv = matrix(nrow = S-1, ncol = Kmax)
                  if("pois" %in% lois){
                    nr = which(lois == "pois")
                    fuv[nr,] = MGM[,1]
                  }
                  if("geom" %in% lois){
                    nr = which(lois == "geom")
                    fuv[nr,] = MGM[,2]
                  }
                  if("nbinom" %in% lois){
                    nr = which(lois == "nbinom")
                    fuv[nr,] = MGM[,3]
                  }
                  if("dweibull" %in% lois){    
                    nr = which(lois == "dweibull")
                    fuv[nr,] = MGM[,4]
                  }
                  if("unif" %in% lois){    
                    nr = which(lois == "unif")
                    fuv[nr,] = MGM[,5]
                  }
                  
                  
                  
                  -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                     # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                     +  sum(vect.Nijk[i+S*(j-1)+plus]*dunif(1:(Kmax), min = 0, max = par[(S*(S-2)+u[1])], log = TRUE))
                     + sum( Nei[i,]* log(1-puv%*%fuv))
                  )  
                  
                  
                }
                
                
                ### CONSTRAINTS
                # puv > 0
                u0 = matrix(0, nrow = (S*(S-2))+2*(S-1), ncol = (S*(S-2))+2*(S-1))
                diag(u0) <- 1
                c0 = rep(0,(S*(S-2))+2*(S-1))
                
                # puv < 1
                u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2))+2*(S-1))
                diag(u1) <- -1
                c1 = c(rep(-1,S*(S-2)))
                
                u2 = rbind(u0,u1)
                c2 = c(c0,c1)
           
                
                ###  Vector for the initialisation of optim
                init = c()
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                  init[p] = c(2, 1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                  init[g] = c(0.5, 0.4)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                  init[nb] = c(2, 4)
                }
                if("dweibull" %in% lois){ 
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                  init[dw] = c(0.6, 0.8)
                }
                if("unif" %in% lois){ 
                  u = 2*which(lois == "unif")-1
                  u = c(u, u+1)
                  init[u] = c(Kmax, Kmax)
                }
                
                
                ### log-likelihood optimisation
                CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), init),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
                ### get the parameters
                parC = CO2$par
                parL = parC[1:(S*(S-2))]
                maxUnif = parC[(S*(S-2)+u[1])]
                
                ### Transformation of parameters
                ## Add parameters of type: pat = 1 - pac - pag
                ## in par vector of optim
                b = 1
                parP = c()
                while( b <= length(parL) ){
                  parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
                  b = b +S-2
                  
                }
                
                ## Add zero in diagonal 
                c = 1 
                parP0 = c()
                while( c < length(parP) ){
                  parP0 = c(parP0, 0,parP[c:(c+(S-1))])
                  c = c + S
                }
                parP0 = c(parP0, 0)
                P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
                # param = array(c(matrix(theta0, nrow = S, ncol = S, byrow = TRUE),matrix(0, nrow = S, ncol = S)), c(S,S,2))
                param[i,j,1] = maxUnif
                
              }
              ## censoring at the beginning
            }else if ( cens.end == 0 && cens.beg == 1 ){
              ########################
              ## Estimation of transition matrix
              #######################
              P = Nij/Ni
              
              plus = seq(from = 0, to = Kmax*S*S-1, by = S*S )
              ## vect.Nijk is the vector of the array Nijk
              vect.Nijk = as.vector(Nijk)
              ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
              ## In the same way, we create the vector  vect.Nbij : 
              ## vect.Nbij is the vector of the array Nbij
              vect.Nbij = as.vector((Nbij))
              ## the Nbij[i,j,] correspond for the vector Nbij to vect.Nbij[i+S*(j-1)+plus]
              ########################
              ## Estimation of parameters
              ########################
              logvrais_theta = function(par){
                -(sum(vect.Nijk[i+S*(j-1)+plus]*dunif(1:(Kmax), min = 0, max = par, log = TRUE))
                  + sum(vect.Nbij[i+S*(j-1)+plus]*log(punif(1:Kmax-1, min = 0, max = par, lower.tail = FALSE)))
                )
              }
              
              out<-optim(par = Kmax, fn = logvrais_theta, method = "Brent", lower = 0, upper = Kmax+1)
              param[i,j,1] = out$par
              
              
              ## censoring at the end and at the beginning
            }else if (cens.end == 1 && cens.beg == 1 ){
              if ( S == 2 ){
                #######################################
                ## Estimation of parameters  ##########
                #######################################
                plus = seq(from = 0, to = Kmax*S*S-1, by = S*S)
                ## vect.Nijk is the vector of the array Nijk
                vect.Nijk = as.vector(Nijk)
                ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
                ## estimation function of  p_{jt+1v} et theta_{jt+1v}
                Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
                El.diag = seq(from= 1, to = S*S, by = S+1)
                Nij.vect = Nij.vect0[-El.diag]
                # Nij0 = Nij0[c(etat.cens),]
                ## vect.Nbij is the vector of the array Nbij
                vect.Nbij = as.vector((Nbij))
                ## the Nbij[i,j,] corresponds for the vector Nbij to vect.Nbij[i+S*(j-1)+plus]
                
                ### Separation of parameters in order to compute log-likelihood
                Nij1 = Nij.vect[1]
                Nij2 = Nij.vect[2]
                lois = distr[i,]
                lois = lois[-i]
                
                ## Function to estimate
                logvrais = function(par){
                  # m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                  # p = cbind(m, 1-rowSums(m))
                  # puv = p[i,]
                  # # ind1 = seq(from = (S*(S-2))+1, to =(S*(S-2)+5), by=2)
                  # ind2 = seq(from = (S*(S-2))+2, to =(S*(S-2)+6), by=2)
                  p = NULL
                  g = NULL
                  nb = NULL
                  dw = NULL
                  u = NULL
                  
                  if("pois" %in% lois){
                    p = 2*which(lois == "pois")-1
                    p = c(p, p+1)
                  }
                  if("geom" %in% lois){
                    g = 2*which(lois == "geom")-1
                    g = c(g, g+1)
                  }
                  if("nbinom" %in% lois){
                    nb = 2*which(lois == "nbinom")-1
                    nb = c(nb, nb+1)
                  }
                  if("dweibull" %in% lois){
                    dw = 2*which(lois == "dweibull")-1
                    dw = c(dw, dw+1)
                  }
                  if("unif" %in% lois){
                    u = 2*which(lois == "unif")-1
                    u = c(u, u+1)
                  }
                  
                  if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                  
                  MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])]), 
                                 dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])]),
                                 dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])]), 
                                 ddw,
                                 dunif(1:Kmax, min = 0, max = par[(S*(S-2)+u[1])])), nrow = Kmax)
                  
                  
                  # print(MGM)
                  fuv = matrix(nrow = S-1, ncol = Kmax)
                  if("pois" %in% lois){
                    nr = which(lois == "pois")
                    fuv[nr,] = MGM[,1]
                  }
                  if("geom" %in% lois){
                    nr = which(lois == "geom")
                    fuv[nr,] = MGM[,2]
                  }
                  if("nbinom" %in% lois){
                    nr = which(lois == "nbinom")
                    fuv[nr,] = MGM[,3]
                  }
                  if("dweibull" %in% lois){    
                    nr = which(lois == "dweibull")
                    fuv[nr,] = MGM[,4]
                  }
                  if("unif" %in% lois){    
                    nr = which(lois == "unif")
                    fuv[nr,] = MGM[,5]
                  }
                  
                  puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                  puv.vect0 = as.vector(t(puv))
                  puv.vect = puv.vect0[-El.diag]
                  
                  -( sum(Nij.vect*log(puv.vect))
                     # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                     +  sum(vect.Nijk[i+S*(j-1)+plus]*dunif(1:(Kmax), min = 0, max = par[(S*(S-2)+u[1])], log = TRUE))
                     + sum(vect.Nbij[i+S*(j-1)+plus]*log(punif(1:Kmax-1, min = 0, max = par[(S*(S-2)+u[1])], lower.tail = FALSE)))
                     + sum( Nei[i,]* log(1-puv.vect%*%fuv))
                  )
                }
                
              
                p = c(0,0)
                g = c(0,0)
                nb = c(0,0)
                dw = c(0,0)
                u = c(0,0)
                
                init = c()
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                  init[p] = c(2, 1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                  init[g] = c(0.5, 0.4)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                  init[nb] = c(2, 4)
                }
                if("dweibull" %in% lois){ 
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                  init[dw] = c(0.6, 0.8)
                }
                if("unif" %in% lois){ 
                  u = 2*which(lois == "unif")-1
                  u = c(u, u+1)
                  init[u] = c(Kmax, Kmax+1)
                }
                
                
                ### Log-likelihood optimisation
                CO2<-optim(par=init,logvrais, method ="Nelder-Mead")
                ### Get the parameterss
                parC = CO2$par
                # parL = parC[1:(S*(S-2))]
                theta = parC[u[1]]
                
                puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                P = puv
                # param = array(c(matrix(theta0, nrow = S, ncol = S, byrow = TRUE),matrix(0, nrow = S, ncol = S)), c(S,S,2))
                param[i,j,1] = theta
                
                
              }else{
                #######################################
                ##Estimation of parameters ##########
                #######################################
                plus = seq(from = 0, to = Kmax*S*S-1, by = S*S )
                ## vect.Nijk is the vector of the array Nijk
                vect.Nijk = as.vector(Nijk)
                ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
                ## In the same way, we create the vector  vect.Nbij : 
                ## vect.Nbij is the vector of the array Nbij
                vect.Nbij = as.vector((Nbij))
                ## the Nbij[i,j,] correspond for the vector Nbij to vect.Nbij[i+S*(j-1)+plus]
                Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
                El.diag = seq(from= 1, to = S*S, by = S+1)
                Nij.vect = Nij.vect0[-El.diag]
                # Nij0 = Nij0[c(etat.cens),]
                
                
                ### Separation of parameters in order to compute the log-likelihood
                Nij1 = c()
                Nij2 = c()
                a = 1
                while (a < length(Nij.vect)){
                  Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## corresponds to the (S-2) first values of Nij on each line
                  Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## corresponds to the last value of Nij on each line
                  a = a + (S-1)
                }
                
                lois = distr[i,]
                lois = lois[-i]
                
               ## Functions to estimate
                logvrais = function(par){
                  m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                  p = cbind(m, 1-rowSums(m))
                  puv = p[i,]
                  
                  p = NULL
                  g = NULL
                  nb = NULL
                  dw = NULL
                  u = NULL
                  
                  # ind1 = seq(from = (S*(S-2))+1, to =(S*(S-2)+5), by=2)
                  # ind2 = seq(from = (S*(S-2))+2, to =(S*(S-2)+6), by=2)
                  if("pois" %in% lois){
                    p = 2*which(lois == "pois")-1
                    p = c(p, p+1)
                  }
                  if("geom" %in% lois){
                    g = 2*which(lois == "geom")-1
                    g = c(g, g+1)
                  }
                  if("nbinom" %in% lois){
                    nb = 2*which(lois == "nbinom")-1
                    nb = c(nb, nb+1)
                  }
                  if("dweibull" %in% lois){
                    dw = 2*which(lois == "dweibull")-1
                    dw = c(dw, dw+1)
                  }
                  if("unif" %in% lois){
                    u = 2*which(lois == "unif")-1
                    u = c(u, u+1)
                  }
                  
                  if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                  
                  MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])]), 
                                 dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])]),
                                 dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])]), 
                                 ddw,
                                 dunif(1:Kmax, min = 0, max = par[S*(S-2)+u[1]])), nrow = Kmax)
                  
                  # print(MGM)
                  fuv = matrix(nrow = S-1, ncol = Kmax)
                  if("pois" %in% lois){
                    nr = which(lois == "pois")
                    fuv[nr,] = MGM[,1]
                  }
                  if("geom" %in% lois){
                    nr = which(lois == "geom")
                    fuv[nr,] = MGM[,2]
                  }
                  if("nbinom" %in% lois){
                    nr = which(lois == "nbinom")
                    fuv[nr,] = MGM[,3]
                  }
                  if("dweibull" %in% lois){    
                    nr = which(lois == "dweibull")
                    fuv[nr,] = MGM[,4]
                  }
                  if("unif" %in% lois){    
                    nr = which(lois == "unif")
                    fuv[nr,] = MGM[,5]
                  }
                  
                  
                  
                  -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                     # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                     +  sum(vect.Nijk[i+S*(j-1)+plus]*dunif(1:(Kmax), min = 0, max = par[(S*(S-2)+u[1])], log = TRUE))
                     + sum(vect.Nbij[i+S*(j-1)+plus]*log(punif(1:Kmax-1, min = 0, max = par[(S*(S-2)+u[1])], lower.tail = FALSE)))
                     + sum( Nei[i,]* log(1-puv%*%fuv))
                  )
                }
                
                
                ### CONSTRAINTS
                # puv > 0
                u0 = matrix(0, nrow = (S*(S-2))+2*(S-1), ncol = (S*(S-2))+2*(S-1))
                diag(u0) <- 1
                c0 = rep(0,(S*(S-2))+2*(S-1))
                
                # puv < 1
                u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2))+2*(S-1))
                diag(u1) <- -1
                c1 = c(rep(-1,S*(S-2)))
                
                u2 = rbind(u0,u1)
                c2 = c(c0,c1)

                p = c(0,0)
                g = c(0,0)
                nb = c(0,0)
                dw = c(0,0)
                u = c(0,0)
                
                init = c()
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                  init[p] = c(2, 1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                  init[g] = c(0.5, 0.4)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                  init[nb] = c(2, 4)
                }
                if("dweibull" %in% lois){ 
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                  init[dw] = c(0.6, 0.8)
                }
                if("unif" %in% lois){ 
                  u = 2*which(lois == "unif")-1
                  u = c(u, u+1)
                  init[u] = c(Kmax, Kmax + 1)
                }
                
                
                ### Log-likelihood optimisation
                CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), init),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
                ### Get the parameters
                parC = CO2$par
                parL = parC[1:(S*(S-2))]
                theta = parC[(S*(S-2)+u[1])]
                
                ### Transformation of parameters
                ## Add parameters of type  : pat = 1 - pac - pag
                ## in par vector of optim
                b = 1
                parP = c()
                while( b <= length(parL) ){
                  parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
                  b = b +S-2
                  
                }
                
               ## Add zeros on diagonal
                c = 1 
                parP0 = c()
                while( c < length(parP) ){
                  parP0 = c(parP0, 0,parP[c:(c+(S-1))])
                  c = c + S
                }
                parP0 = c(parP0, 0)
                P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
                # param = array(c(matrix(theta0, nrow = S, ncol = S, byrow = TRUE),matrix(0, nrow = S, ncol = S)), c(S,S,2))
                param[i,j,1] = theta
                
              }
              
              
            }
          }else if ( distr[i,j] == "pois" ){
            ## censoring at the end
            if( cens.end == 1 && cens.beg == 0){
              
              ## Particular case S = 2 
               if ( S == 2 ){
                #############################
                ## Estimation of parameters 
                #############################
                plus = seq(from = 0, to = Kmax*S*S-1, by = S*S)
                ## vect.Nijk is the vector of the array Nijk
                vect.Nijk = as.vector(Nijk)
                ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
                ## estimation function of  p_{jt+1v} et theta_{jt+1v}
                Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
                El.diag = seq(from= 1, to = S*S, by = S+1)
                Nij.vect = Nij.vect0[-El.diag]
                # Nij0 = Nij0[c(etat.cens),]
                
                
                ### Separation of parameters in order to compute log-likelihood
                Nij1 = Nij.vect[1]
                Nij2 = Nij.vect[2]
                lois = distr[i,]
                lois = lois[-i]
                
               ## Functions to estimate
                logvrais = function(par){
                  p = NULL
                  g = NULL
                  nb = NULL
                  dw = NULL
                  
                  if("pois" %in% lois){
                    p = 2*which(lois == "pois")-1
                    p = c(p, p+1)
                  }
                  if("geom" %in% lois){
                    g = 2*which(lois == "geom")-1
                    g = c(g, g+1)
                  }
                  if("nbinom" %in% lois){
                    nb = 2*which(lois == "nbinom")-1
                    nb = c(nb, nb+1)
                  }
                  if("dweibull" %in% lois){
                    dw = 2*which(lois == "dweibull")-1
                    dw = c(dw, dw+1)
                  }
                  
                  if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                  
                  MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])]), 
                                 dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])]),
                                 dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])]), 
                                 ddw), nrow = Kmax)
                  
                  
                  # print(MGM)
                  fuv = matrix(nrow = S-1, ncol = Kmax)
                  if("pois" %in% lois){
                    nr = which(lois == "pois")
                    fuv[nr,] = MGM[,1]
                  }
                  if("geom" %in% lois){
                    nr = which(lois == "geom")
                    fuv[nr,] = MGM[,2]
                  }
                  if("nbinom" %in% lois){
                    nr = which(lois == "nbinom")
                    fuv[nr,] = MGM[,3]
                  }
                  if("dweibull" %in% lois){    
                    nr = which(lois == "dweibull")
                    fuv[nr,] = MGM[,4]
                  }
                  
                  puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                  puv.vect0 = as.vector(t(puv))
                  puv.vect = puv.vect0[-El.diag]
                  
                  -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                     # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                     +  sum(vect.Nijk[i+S*(j-1)+plus]*dpois(0:(Kmax-1), lambda = par[p[1]], log = TRUE))
                     + sum( Nei[i,]* log(1-puv.vect%*%fuv))
                  )  
                  
                  
                }
                
                
                ### Vector to initialise in optim
                init = c()
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                  init[p] = c(2, 1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                  init[g] = c(0.5, 0.4)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                  init[nb] = c(2, 4)
                }
                if("dweibull" %in% lois){ 
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                  init[dw] = c(0.6, 0.8)
                }
                
                
                ### Log-likelihood optimisation
                CO2<-optim(par = init,fn = logvrais, method ="Nelder-Mead")
                ### Get the parameters
                parC = CO2$par
                theta = parC[p[1]]
                param[i,j,1] = theta
                
                puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                P = puv
                
              }else{
                #############################
                ## Estimation of parameters 
                #############################
                plus = seq(from = 0, to = Kmax*S*S-1, by = S*S)
                ## vect.Nijk is the vector of the array Nijk
                vect.Nijk = as.vector(Nijk)
                ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
                ## estimation function of  p_{jt+1v} et theta_{jt+1v}
                Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
                El.diag = seq(from= 1, to = S*S, by = S+1)
                Nij.vect = Nij.vect0[-El.diag]
                # Nij0 = Nij0[c(etat.cens),]
                
                
                ### Separation of parameters in order to compute the log-likelihood
                Nij1 = c()
                Nij2 = c()
                a = 1
                while (a < length(Nij.vect)){
                  Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
                  Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
                  a = a + (S-1)
                }
                
                lois = distr[i,]
                lois = lois[-i]
                
               ## Functions to estimate
                logvrais = function(par){
                  m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                  p = cbind(m, 1-rowSums(m))
                  puv = p[i,]
                  # ind1 = seq(from = (S*(S-2))+1, to =(S*(S-2)+5), by=2)
                  # ind2 = seq(from = (S*(S-2))+2, to =(S*(S-2)+6), by=2)
                  p = NULL
                  g = NULL
                  nb = NULL
                  dw = NULL
                  
                  if("pois" %in% lois){
                    p = 2*which(lois == "pois")-1
                    p = c(p, p+1)
                  }
                  if("geom" %in% lois){
                    g = 2*which(lois == "geom")-1
                    g = c(g, g+1)
                  }
                  if("nbinom" %in% lois){
                    nb = 2*which(lois == "nbinom")-1
                    nb = c(nb, nb+1)
                  }
                  if("dweibull" %in% lois){
                    dw = 2*which(lois == "dweibull")-1
                    dw = c(dw, dw+1)
                  }
                  
                  if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                  
                  MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])]), 
                                 dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])]),
                                 dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])]), 
                                 ddw), nrow = Kmax)
                  
                  # print(MGM)
                  fuv = matrix(nrow = S-1, ncol = Kmax)
                  if("pois" %in% lois){
                    nr = which(lois == "pois")
                    fuv[nr,] = MGM[,1]
                  }
                  if("geom" %in% lois){
                    nr = which(lois == "geom")
                    fuv[nr,] = MGM[,2]
                  }
                  if("nbinom" %in% lois){
                    nr = which(lois == "nbinom")
                    fuv[nr,] = MGM[,3]
                  }
                  if("dweibull" %in% lois){    
                    nr = which(lois == "dweibull")
                    fuv[nr,] = MGM[,4]
                  }
                  
                  
                  
                  -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                     # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                     +  sum(vect.Nijk[i+S*(j-1)+plus]*dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])], log = TRUE))
                     + sum( Nei[i,]* log(1-puv%*%fuv))
                  )  
                  
                  
                }
                
                
                ### CONSTRAINTS
                # puv > 0
                u0 = matrix(0, nrow = (S*(S-2))+2*(S-1), ncol = (S*(S-2))+2*(S-1))
                diag(u0) <- 1
                c0 = rep(0,(S*(S-2))+2*(S-1))
                
                # puv < 1
                u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2))+2*(S-1))
                diag(u1) <- -1
                c1 = c(rep(-1,S*(S-2)))
                
                u2 = rbind(u0,u1)
                c2 = c(c0,c1)
                
            
                ### Vector to initialise par in optim
                init = c()
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                  init[p] = c(2, 1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                  init[g] = c(0.5, 0.4)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                  init[nb] = c(2, 4)
                }
                if("dweibull" %in% lois){ 
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                  init[dw] = c(0.6, 0.8)
                }
                
                
                ### Log-likelihood optimisation
                CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), init),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
                ### Get the parameters
                parC = CO2$par
                parL = parC[1:(S*(S-2))]
                theta = parC[(S*(S-2)+p[1])]
                              
                ### Transformation of parameters
                ## Add parameters of type: pat = 1 - pac - pag
                ## in par vector of optim
                b = 1
                parP = c()
                while( b <= length(parL) ){
                  parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
                  b = b +S-2
                  
                }
                
               ## Add zeros in diagonale 
                c = 1 
                parP0 = c()
                while( c < length(parP) ){
                  parP0 = c(parP0, 0,parP[c:(c+(S-1))])
                  c = c + S
                }
                parP0 = c(parP0, 0)
                P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
                # param = array(c(matrix(theta0, nrow = S, ncol = S, byrow = TRUE),matrix(0, nrow = S, ncol = S)), c(S,S,2))
                param[i,j,1] = theta
            
              }    
                ##  censoring at the beginning 
            }else if( cens.end == 0 && cens.beg == 1 ){   
              
              ########################
              ## Estimation of the transition matrix
              #######################
              P = Nij/Ni
              
              plus = seq(from = 0, to = Kmax*S*S-1, by = S*S )
              ## vect.Nijk is the vector of the array Nijk
              vect.Nijk = as.vector(Nijk)
              ##  Nijk[i,j,] corresponds for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
              ## In the same way, we create the vector  vect.Nbij : 
              ## vect.Nbij is the vector of the array Nbij
              vect.Nbij = as.vector((Nbij))
              ## the Nbij[i,j,] correspond for the vector Nbij to vect.Nbij[i+S*(j-1)+plus]
              ########################
              ## Estimation of parameters  of the distribution 
              ########################
              logvrais_theta = function(par){
                -(sum(vect.Nijk[i+S*(j-1)+plus]*dpois(0:(Kmax-1), lambda = par, log = TRUE))
                  + sum(vect.Nbij[i+S*(j-1)+plus]*log(ppois(1:Kmax-1, lambda = par, lower.tail = FALSE)))
                )
              }
              
              out<-optim(par = 1, fn = logvrais_theta, method = "Brent", lower = 0, upper = Kmax)
              param[i,j,1] = out$par
              
              ## censoring at the end and at the beginning
            }else if( cens.end == 1 && cens.beg == 1 ){
              
              if ( S == 2 ){
                #######################################
                ##Estimation of parameters ##########
                #######################################
                plus = seq(from = 0, to = Kmax*S*S-1, by = S*S)
                ## vect.Nijk is the vector of the array Nijk
                vect.Nijk = as.vector(Nijk)
                ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
                ## estimation function of  p_{jt+1v} and theta_{jt+1v}
                Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
                El.diag = seq(from= 1, to = S*S, by = S+1)
                Nij.vect = Nij.vect0[-El.diag]
                # Nij0 = Nij0[c(etat.cens),]
                ## vect.Nbij is the vector of the array Nbij
                vect.Nbij = as.vector((Nbij))
                ## the Nbij[i,j,] correspond for the vector Nbij to vect.Nbij[i+S*(j-1)+plus]
                
                
                ### Separation of parameters in order to compute the log-likelihood
                Nij1 = Nij.vect[1]
                Nij2 = Nij.vect[2]
                # a = 1
                # while (a < length(Nij.vect)){
                #   Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
                #   Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
                #   a = a + (S-1)
                # }
                # 
                lois = distr[i,]
                lois = lois[-i]
                
               ## Functions to estimate
                logvrais = function(par){
                  # m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                  # p = cbind(m, 1-rowSums(m))
                  # puv = p[i,]
                  # # ind1 = seq(from = (S*(S-2))+1, to =(S*(S-2)+5), by=2)
                  # ind2 = seq(from = (S*(S-2))+2, to =(S*(S-2)+6), by=2)
                  p = NULL
                  g = NULL
                  nb = NULL
                  dw = NULL
                  
                  if("pois" %in% lois){
                    p = 2*which(lois == "pois")-1
                    p = c(p, p+1)
                  }
                  if("geom" %in% lois){
                    g = 2*which(lois == "geom")-1
                    g = c(g, g+1)
                  }
                  if("nbinom" %in% lois){
                    nb = 2*which(lois == "nbinom")-1
                    nb = c(nb, nb+1)
                  }
                  if("dweibull" %in% lois){
                    dw = 2*which(lois == "dweibull")-1
                    dw = c(dw, dw+1)
                  }
                  
                  if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                  
                  MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])]), 
                                 dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])]),
                                 dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])]), 
                                 ddw), nrow = Kmax)
                  
                  
                  # print(MGM)
                  fuv = matrix(nrow = S-1, ncol = Kmax)
                  if("pois" %in% lois){
                    nr = which(lois == "pois")
                    fuv[nr,] = MGM[,1]
                  }
                  if("geom" %in% lois){
                    nr = which(lois == "geom")
                    fuv[nr,] = MGM[,2]
                  }
                  if("nbinom" %in% lois){
                    nr = which(lois == "nbinom")
                    fuv[nr,] = MGM[,3]
                  }
                  if("dweibull" %in% lois){    
                    nr = which(lois == "dweibull")
                    fuv[nr,] = MGM[,4]
                  }
                  
                  puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                  puv.vect0 = as.vector(t(puv))
                  puv.vect = puv.vect0[-El.diag]
                  
                  -( sum(Nij.vect*log(puv.vect))
                     # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                     +  sum(vect.Nijk[i+S*(j-1)+plus]*dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])], log = TRUE))
                     + sum(vect.Nbij[i+S*(j-1)+plus]*log(ppois(1:Kmax-1, lambda = par[(S*(S-2)+p[1])], lower.tail = FALSE)))
                     + sum( Nei[i,]* log(1-puv.vect%*%fuv))
                  )
                }
                
                
                ### CONSTRAINTS
                # puv > 0
                # u0 = matrix(0, nrow = (S*(S-2))+2*(S-1), ncol = (S*(S-2))+2*(S-1))
                # diag(u0) <- 1
                # c0 = rep(0,(S*(S-2))+2*(S-1))
                # 
                # # puv < 1
                # u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2))+2*(S-1))
                # diag(u1) <- -1
                # c1 = c(rep(-1,S*(S-2)))
                # 
                # u2 = rbind(u0,u1)
                # c2 = c(c0,c1)
                # 
                ## verification of constraints
                # k2 <- length(c2)
                # n <- dim(u2)[2]
                # for(i in seq_len(k2)) {
                #   j <- which( u2[i,] != 0 )
                #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
                #   cat(" >= " )
                #   cat( c2[i], "\n" )
                # }
                p = c(0,0)
                g = c(0,0)
                nb = c(0,0)
                dw = c(0,0)
                
                init = c()
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                  init[p] = c(2, 1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                  init[g] = c(0.5, 0.4)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                  init[nb] = c(2, 4)
                }
                if("dweibull" %in% lois){ 
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                  init[dw] = c(0.6, 0.8)
                }
                
                
                ### Log-likelihood optimisation
                CO2<-optim(par=init,logvrais, method ="Nelder-Mead")
                ### Get the parameters
                parC = CO2$par
                # parL = parC[1:(S*(S-2))]
                theta = parC[p[1]]
                
                puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                P = puv
                # param = array(c(matrix(theta0, nrow = S, ncol = S, byrow = TRUE),matrix(0, nrow = S, ncol = S)), c(S,S,2))
                param[i,j,1] = theta
                
                
              }else{
                #######################################
                ##Estimation of parameters ##########
                #######################################
                plus = seq(from = 0, to = Kmax*S*S-1, by = S*S )
                ## vect.Nijk is the vector of the array Nijk
                vect.Nijk = as.vector(Nijk)
                ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
                ## In the same way, we create the vector  vect.Nbij : 
                ## vect.Nbij is the vector of the array Nbij
                vect.Nbij = as.vector((Nbij))
                ##  Nbij[i,j,] corresponds for the vector Nbij to vect.Nbij[i+S*(j-1)+plus]
                Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
                El.diag = seq(from= 1, to = S*S, by = S+1)
                Nij.vect = Nij.vect0[-El.diag]
                # Nij0 = Nij0[c(etat.cens),]
                
                
                ### Separation of parameters in order to compute log-likelihood
                Nij1 = c()
                Nij2 = c()
                a = 1
                while (a < length(Nij.vect)){
                  Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
                  Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
                  a = a + (S-1)
                }
                
                lois = distr[i,]
                lois = lois[-i]
                
               ## Functions to estimate
                logvrais = function(par){
                  m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                  p = cbind(m, 1-rowSums(m))
                  puv = p[i,]
                  
                  p = NULL
                  g = NULL
                  nb = NULL
                  dw = NULL
                  
                  # ind1 = seq(from = (S*(S-2))+1, to =(S*(S-2)+5), by=2)
                  # ind2 = seq(from = (S*(S-2))+2, to =(S*(S-2)+6), by=2)
                  if("pois" %in% lois){
                    p = 2*which(lois == "pois")-1
                    p = c(p, p+1)
                  }
                  if("geom" %in% lois){
                    g = 2*which(lois == "geom")-1
                    g = c(g, g+1)
                  }
                  if("nbinom" %in% lois){
                    nb = 2*which(lois == "nbinom")-1
                    nb = c(nb, nb+1)
                  }
                  if("dweibull" %in% lois){
                    dw = 2*which(lois == "dweibull")-1
                    dw = c(dw, dw+1)
                  }
                  if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                  
                  MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])]), 
                                 dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])]),
                                 dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])]), 
                                 ddw), nrow = Kmax)
                  # print(MGM)
                  fuv = matrix(nrow = S-1, ncol = Kmax)
                  if("pois" %in% lois){
                    nr = which(lois == "pois")
                    fuv[nr,] = MGM[,1]
                  }
                  if("geom" %in% lois){
                    nr = which(lois == "geom")
                    fuv[nr,] = MGM[,2]
                  }
                  if("nbinom" %in% lois){
                    nr = which(lois == "nbinom")
                    fuv[nr,] = MGM[,3]
                  }
                  if("dweibull" %in% lois){    
                    nr = which(lois == "dweibull")
                    fuv[nr,] = MGM[,4]
                  }
                  
                  
                  
                  -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                     # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                     +  sum(vect.Nijk[i+S*(j-1)+plus]*dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])], log = TRUE))
                     + sum(vect.Nbij[i+S*(j-1)+plus]*log(ppois(1:Kmax-1, lambda = par[(S*(S-2)+p[1])], lower.tail = FALSE)))
                     + sum( Nei[i,]* log(1-puv%*%fuv))
                  )
                }
                
                
                ### CONSTRAINTS
                # puv > 0
                u0 = matrix(0, nrow = (S*(S-2))+2*(S-1), ncol = (S*(S-2))+2*(S-1))
                diag(u0) <- 1
                c0 = rep(0,(S*(S-2))+2*(S-1))
                
                # puv < 1
                u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2))+2*(S-1))
                diag(u1) <- -1
                c1 = c(rep(-1,S*(S-2)))
                
                u2 = rbind(u0,u1)
                c2 = c(c0,c1)
                
                ## verification of constraints
                # k2 <- length(c2)
                # n <- dim(u2)[2]
                # for(i in seq_len(k2)) {
                #   j <- which( u2[i,] != 0 )
                #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
                #   cat(" >= " )
                #   cat( c2[i], "\n" )
                # }
                p = c(0,0)
                g = c(0,0)
                nb = c(0,0)
                dw = c(0,0)
                
                init = c()
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                  init[p] = c(2, 1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                  init[g] = c(0.5, 0.4)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                  init[nb] = c(2, 4)
                }
                if("dweibull" %in% lois){ 
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                  init[dw] = c(0.6, 0.8)
                }
                
                
                ### Log-likelihood optimisation
                CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), init),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
                ### Get the parameters
                parC = CO2$par
                parL = parC[1:(S*(S-2))]
                theta = parC[(S*(S-2)+p[1])]
                
                ### Transformation of parameters
                ## Add parameters of type  : pat = 1 - pac - pag
                ## in par vector of optim
                b = 1
                parP = c()
                while( b <= length(parL) ){
                  parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
                  b = b +S-2
                  
                }
                
               ## Add zeros in diagonale 
                c = 1 
                parP0 = c()
                while( c < length(parP) ){
                  parP0 = c(parP0, 0,parP[c:(c+(S-1))])
                  c = c + S
                }
                parP0 = c(parP0, 0)
                P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
                # param = array(c(matrix(theta0, nrow = S, ncol = S, byrow = TRUE),matrix(0, nrow = S, ncol = S)), c(S,S,2))
                param[i,j,1] = theta
                
              }
                            
              ## No censoring  
            }else{        
              
              ##############
              ## Estimation of transition matrix
              ##############
              P = Nij/Ni
              P[which(is.na(P))] = 0
              
              ##############
              ## Estimation of theta : 
              #############
              ## vector plus for getting all the Nij for each k
              plus = seq(from = 0, to = Kmax*S*S-1, by = S*S )
              ## vect.Nijk is the vector of the array Nijk
              vect.Nijk = as.vector(Nijk)
              ##  Nijk[i,j,] corresponds for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
              logvrais_theta = function(par){
                -(sum(vect.Nijk[i+S*(j-1)+plus]*dpois(0:(Kmax-1), lambda = par, log = TRUE)))
              }
              
              if (i != j ){ 
                ## Warning for warnings due to par initialization
                # out<-optim(par = 1, logvrais_theta, method = "BFGS")
                # when maximization is posssible withot optim, do without
                x = rep(0:(Kmax-1), vect.Nijk[i+S*(j-1)+plus])
                param[i,j,1] = (sum(x)/sum(vect.Nijk[i+S*(j-1)+plus]))
                # param[i,j,1] = out$par
              }
              
            } 
            
            
            
          }else if ( distr[i,j] == "geom" ){
            
            ## censoring at the end
            if( cens.end == 1 && cens.beg == 0){
              if (S == 2){
                #############################
                ## Estimation of parameters 
                #############################
                plus = seq(from = 0, to = Kmax*S*S-1, by = S*S)
                ## vect.Nijk is the vector of the array Nijk
                vect.Nijk = as.vector(Nijk)
                ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
                ## estimation function of  p_{jt+1v} and theta_{jt+1v}
                Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
                El.diag = seq(from= 1, to = S*S, by = S+1)
                Nij.vect = Nij.vect0[-El.diag]
                # Nij0 = Nij0[c(etat.cens),]
                
                
                ### Separation of parameters in order to compute the log-likelihood
                Nij1 = Nij.vect[1]
                Nij2 = Nij.vect[2]
                # a = 1
                # while (a < length(Nij.vect)){
                #   Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## corresponds to the (S-2) first values of Nij on each line
                #   Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## corresponds to the last value of Nij on each line
                #   a = a + (S-1)
                # }
                # 
                lois = distr[i,]
                lois = lois[-i]
                
               ## Functions to estimate
                logvrais = function(par){
                  # m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                  # p = cbind(m, 1-rowSums(m))
                  # puv = p[i,]
                  # # ind1 = seq(from = (S*(S-2))+1, to =(S*(S-2)+5), by=2)
                  # ind2 = seq(from = (S*(S-2))+2, to =(S*(S-2)+6), by=2)
                  p = NULL
                  g = NULL
                  nb = NULL
                  dw = NULL
                  
                  if("pois" %in% lois){
                    p = 2*which(lois == "pois")-1
                    p = c(p, p+1)
                  }
                  if("geom" %in% lois){
                    g = 2*which(lois == "geom")-1
                    g = c(g, g+1)
                  }
                  if("nbinom" %in% lois){
                    nb = 2*which(lois == "nbinom")-1
                    nb = c(nb, nb+1)
                  }
                  if("dweibull" %in% lois){
                    dw = 2*which(lois == "dweibull")-1
                    dw = c(dw, dw+1)
                  }
                  
                  if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                  
                  MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])]), 
                                 dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])]),
                                 dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])]), 
                                 ddw), nrow = Kmax)
                  
                  
                  # print(MGM)
                  fuv = matrix(nrow = S-1, ncol = Kmax)
                  if("pois" %in% lois){
                    nr = which(lois == "pois")
                    fuv[nr,] = MGM[,1]
                  }
                  if("geom" %in% lois){
                    nr = which(lois == "geom")
                    fuv[nr,] = MGM[,2]
                  }
                  if("nbinom" %in% lois){
                    nr = which(lois == "nbinom")
                    fuv[nr,] = MGM[,3]
                  }
                  if("dweibull" %in% lois){    
                    nr = which(lois == "dweibull")
                    fuv[nr,] = MGM[,4]
                  }
                  
                  puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                  puv.vect0 = as.vector(t(puv))
                  puv.vect = puv.vect0[-El.diag]
                  
                  -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                     # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                     +  sum(vect.Nijk[i+S*(j-1)+plus]*dgeom(0:(Kmax-1), prob = par[g[1]], log = TRUE))
                     + sum( Nei[i,]* log(1-puv.vect%*%fuv))
                  )  
                  
                  
                }
                
                
                ### CONSTRAINTS
                # # puv > 0
                # u0 = matrix(0, nrow = (S*(S-2))+2*(S-1), ncol = (S*(S-2))+2*(S-1))
                # diag(u0) <- 1
                # c0 = rep(0,(S*(S-2))+2*(S-1))
                # 
                # # puv < 1
                # u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2))+2*(S-1))
                # diag(u1) <- -1
                # c1 = c(rep(-1,S*(S-2)))
                # 
                # u2 = rbind(u0,u1)
                # c2 = c(c0,c1)
                
                ## verification of constraints
                # k2 <- length(c2)
                # n <- dim(u2)[2]
                # for(i in seq_len(k2)) {
                #   j <- which( u2[i,] != 0 )
                #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
                #   cat(" >= " )
                #   cat( c2[i], "\n" )
                # }
                
                
                ### Vector to initialise par in optim
                init = c()
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                  init[p] = c(2, 1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                  init[g] = c(0.5, 0.4)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                  init[nb] = c(2, 4)
                }
                if("dweibull" %in% lois){ 
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                  init[dw] = c(0.6, 0.8)
                }
                
                
                ### Log-likelihood optimisation
                CO2<-optim(par = init,fn = logvrais, method ="Nelder-Mead")
                ### Get the parameters
                parC = CO2$par
                theta = parC[g[1]]
                
                puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                P = puv
                param[i,j,1] = theta
                
              }else{
                #############################
                ## Estimation of parameters 
                #############################
                # tNijk = aperm(Nijk, c(2,1,3))i
                # Nijk.vect = as.vector(tNijk)
                # n = 1
                # d = 1
                # Diag = c()
                # while ( n < S*S*Kmax ){
                #   Diag = c(Diag, seq(from = n, to = d*S*S, by = S+1))
                #   d=d+1
                #   n=n+S*S
                # }
                # Nijk.vect = Nijk.vect[-Diag]
                plus = seq(from = 0, to = Kmax*S*S-1, by = S*S)
                ## vect.Nijk is the vector of the array Nijk
                vect.Nijk = as.vector(Nijk)
                ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
                ## estimation function of  p_{jt+1v} and theta_{jt+1v}
                Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
                El.diag = seq(from= 1, to = S*S, by = S+1)
                Nij.vect = Nij.vect0[-El.diag]
                # Nij0 = Nij0[c(etat.cens),]
                
                # 
                ### Separation of parameters in order to compute log-likelihood
                Nij1 = c()
                Nij2 = c()
                a = 1
                while (a < length(Nij.vect)){
                  Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## corresponds to the (S-2) first values of Nij on each line
                  Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## corresponds to the last value of Nij on each line
                  a = a + (S-1)
                }
                
                lois = distr[i,]
                lois = lois[-i]
                
               ## Functions to estimate
                logvrais = function(par){
                  m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                  p = cbind(m, 1-rowSums(m))
                  puv = p[i,]
                  p = NULL
                  g = NULL
                  nb = NULL
                  dw = NULL
                  
                  # ind1 = seq(from = (S*(S-2))+1, to =(S*(S-2)+5), by=2)
                  # ind2 = seq(from = (S*(S-2))+2, to =(S*(S-2)+6), by=2)
                  if("pois" %in% lois){
                    p = 2*which(lois == "pois")-1
                    p = c(p, p+1)
                  }
                  if("geom" %in% lois){
                    g = 2*which(lois == "geom")-1
                    g = c(g, g+1)
                  }
                  if("nbinom" %in% lois){
                    nb = 2*which(lois == "nbinom")-1
                    nb = c(nb, nb+1)
                  }
                  if("dweibull" %in% lois){
                    dw = 2*which(lois == "dweibull")-1
                    dw = c(dw, dw+1)
                  }
                  
                  if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                  
                  MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])]), 
                                 dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])]),
                                 dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])]), 
                                 ddw), nrow = Kmax)
                  
                  # print(MGM)
                  fuv = matrix(nrow = S-1, ncol = Kmax)
                  if("pois" %in% lois){
                    nr = which(lois == "pois")
                    fuv[nr,] = MGM[,1]
                  }
                  if("geom" %in% lois){
                    nr = which(lois == "geom")
                    fuv[nr,] = MGM[,2]
                  }
                  if("nbinom" %in% lois){
                    nr = which(lois == "nbinom")
                    fuv[nr,] = MGM[,3]
                  }
                  if("dweibull" %in% lois){    
                    nr = which(lois == "dweibull")
                    fuv[nr,] = MGM[,4]
                  }
                  -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                     # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                     +  sum(vect.Nijk[i+S*(j-1)+plus]*dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])], log = TRUE))
                     + sum( Nei[i,]* log(1-puv%*%fuv))
                  )
                }
                
                
                ### CONSTRAINTS
                # puv > 0
                u0 = matrix(0, nrow = (S*(S-2))+2*(S-1), ncol = (S*(S-2))+2*(S-1))
                diag(u0) <- 1
                c0 = rep(0,(S*(S-2))+2*(S-1))
                
                # puv < 1
                u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2))+2*(S-1))
                diag(u1) <- -1
                c1 = c(rep(-1,S*(S-2)))
                
                u2 = rbind(u0,u1)
                c2 = c(c0,c1)
                
                ## verification of constraints
                # k2 <- length(c2)
                # n <- dim(u2)[2]
                # for(i in seq_len(k2)) {
                #   j <- which( u2[i,] != 0 )
                #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
                #   cat(" >= " )
                #   cat( c2[i], "\n" )
                # }
                p = c(0,0)
                g = c(0,0)
                nb = c(0,0)
                dw = c(0,0)
                
                init = c()
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                  init[p] = c(2, 1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                  init[g] = c(0.5, 0.4)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                  init[nb] = c(2, 4)
                }
                if("dweibull" %in% lois){ 
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                  init[dw] = c(0.6, 0.8)
                }
                
                ### Log-likelihood optimisation
                CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), init),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
                ### Get the parameters
                parC = CO2$par
                parL = parC[1:(S*(S-2))]
                theta = parC[(S*(S-2)+g[1])]
                
                ### Transformation of parameters
                ## Add parameters of type: pat = 1 - pac - pag
                ## in par vector of optim
                b = 1
                parP = c()
                while( b <= length(parL) ){
                  parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
                  b = b +S-2
                  
                }
                
               ## Add zeros in diagonale 
                c = 1 
                parP0 = c()
                while( c < length(parP) ){
                  parP0 = c(parP0, 0,parP[c:(c+(S-1))])
                  c = c + S
                }
                parP0 = c(parP0, 0)
                P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
                P[which(is.na(P))] = 0
                # param = array(c(matrix(theta0, nrow = S, ncol = S, byrow = TRUE),matrix(0, nrow = S, ncol = S)), c(S,S,2))
                param[i,j,1] = theta
              }
            ##  censoring at the beginning 
          }else if( cens.end == 0 && cens.beg == 1 ){   
            
            ########################
            ## Estimation of the transition matrix
            #######################
            P = Nij/Ni
            P[which(is.na(P))] = 0
            
            plus = seq(from = 0, to = Kmax*S*S-1, by = S*S )
            ## vect.Nijk is the vector of the array Nijk
            vect.Nijk = as.vector(Nijk)
            ##  Nijk[i,j,] correspondS for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
            ## In the same way, we create the vector  vect.Nbij : 
            ## vect.Nbij is the vector of the array Nbij
            vect.Nbij = as.vector((Nbij))
            ## the Nbij[i,j,] corresponds for the vector Nbij to vect.Nbij[i+S*(j-1)+plus]
            ########################
            ## Estimation of parameters  of the distribution 
            ########################
            logvrais_theta = function(par){
              -(sum(vect.Nijk[i+S*(j-1)+plus]*dgeom(0:(Kmax-1), prob = par, log = TRUE))
                + sum(vect.Nbij[i+S*(j-1)+plus]*log(pgeom(1:Kmax-1, prob = par, lower.tail = FALSE)))
              )
            }
            
            out<-optim(par = 0.5, logvrais_theta, method = "Brent", lower = 0, upper = 1)
            param[i,j,1] = out$par
            
            ## censoring at the end and at the beginning
          }else if( cens.end == 1 && cens.beg == 1 ){
            if ( S == 2 ){
              #######################################
              ##Estimation of parameters ##########
              #######################################
              plus = seq(from = 0, to = Kmax*S*S-1, by = S*S)
              ## vect.Nijk is the vector of the array Nijk
              vect.Nijk = as.vector(Nijk)
              ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
              ## estimation function of  p_{jt+1v} and theta_{jt+1v}
              Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
              El.diag = seq(from= 1, to = S*S, by = S+1)
              Nij.vect = Nij.vect0[-El.diag]
              # Nij0 = Nij0[c(etat.cens),]
              ## vect.Nbij is the vector of the array Nbij
              vect.Nbij = as.vector((Nbij))
              ## the Nbij[i,j,] correspond for the vector Nbij to vect.Nbij[i+S*(j-1)+plus]
              
              
              ### Separation of parameters for log-likelihood estimation
              Nij1 = Nij.vect[1]
              Nij2 = Nij.vect[2]
              # a = 1
              # while (a < length(Nij.vect)){
              #   Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## corresponds to the (S-2) first values of Nij on each line
              #   Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## corresponds to the last value of Nij on each line
              #   a = a + (S-1)
              # }
              # 
              lois = distr[i,]
              lois = lois[-i]
              
             ## Functions to estimate
              logvrais = function(par){
                # m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                # p = cbind(m, 1-rowSums(m))
                # puv = p[i,]
                # # ind1 = seq(from = (S*(S-2))+1, to =(S*(S-2)+5), by=2)
                # ind2 = seq(from = (S*(S-2))+2, to =(S*(S-2)+6), by=2)
                p = NULL
                g = NULL
                nb = NULL
                dw = NULL
                
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                }
                if("dweibull" %in% lois){
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                }
                
                if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                
                MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])]), 
                               dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])]),
                               dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])]), 
                               ddw), nrow = Kmax)
                
                
                # print(MGM)
                fuv = matrix(nrow = S-1, ncol = Kmax)
                if("pois" %in% lois){
                  nr = which(lois == "pois")
                  fuv[nr,] = MGM[,1]
                }
                if("geom" %in% lois){
                  nr = which(lois == "geom")
                  fuv[nr,] = MGM[,2]
                }
                if("nbinom" %in% lois){
                  nr = which(lois == "nbinom")
                  fuv[nr,] = MGM[,3]
                }
                if("dweibull" %in% lois){    
                  nr = which(lois == "dweibull")
                  fuv[nr,] = MGM[,4]
                }
                
                puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                puv.vect0 = as.vector(t(puv))
                puv.vect = puv.vect0[-El.diag]
                
                -( sum(Nij.vect*log(puv.vect))
                   # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                   +  sum(vect.Nijk[i+S*(j-1)+plus]*dgeom(0:(Kmax-1), prob = par[g[1]], log = TRUE))
                   + sum(vect.Nbij[i+S*(j-1)+plus]*log(pgeom(1:Kmax-1, prob = par[g[1]], lower.tail = FALSE)))
                   + sum( Nei[i,]* log(1-puv.vect%*%fuv))
                )
              }
              
              
              ### CONSTRAINTS
               p = c(0,0)
              g = c(0,0)
              nb = c(0,0)
              dw = c(0,0)
              
              init = c()
              if("pois" %in% lois){
                p = 2*which(lois == "pois")-1
                p = c(p, p+1)
                init[p] = c(2, 1)
              }
              if("geom" %in% lois){
                g = 2*which(lois == "geom")-1
                g = c(g, g+1)
                init[g] = c(0.5, 0.4)
              }
              if("nbinom" %in% lois){
                nb = 2*which(lois == "nbinom")-1
                nb = c(nb, nb+1)
                init[nb] = c(2, 4)
              }
              if("dweibull" %in% lois){ 
                dw = 2*which(lois == "dweibull")-1
                dw = c(dw, dw+1)
                init[dw] = c(0.6, 0.8)
              }
              
              
              ### Log-likelihood optimisation
              CO2<-optim(par=init,logvrais, method ="Nelder-Mead")
              ### Get the parameters
              parC = CO2$par
              # parL = parC[1:(S*(S-2))]
              theta = parC[g[1]]
              
              puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
              P = puv
              param[i,j,1] = theta
              
              
            }else{
              #############################
              ## Estimation of parameters 
              #############################
              plus = seq(from = 0, to = Kmax*S*S-1, by = S*S )
              ## vect.Nijk is the vector of the array Nijk
              vect.Nijk = as.vector(Nijk)
              ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
              ## In the same way, we create the vector  vect.Nbij : 
              ## vect.Nbij is the vector of the array Nbij
              vect.Nbij = as.vector((Nbij))
              ## the Nbij[i,j,] correspond for the vector Nbij to vect.Nbij[i+S*(j-1)+plus]
              Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
              El.diag = seq(from= 1, to = S*S, by = S+1)
              Nij.vect = Nij.vect0[-El.diag]
              # Nij0 = Nij0[c(etat.cens),]

                
              ### Separation of parameters for log-likelihood estimation
              Nij1 = c()
              Nij2 = c()
              a = 1
              while (a < length(Nij.vect)){
                Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
                Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
                a = a + (S-1)
              }
              
              lois = distr[i,]
              lois = lois[-i]
              
             ## Functions to estimate
              logvrais = function(par){
                m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                p = cbind(m, 1-rowSums(m))
                puv = p[i,]
                p = NULL
                g = NULL
                nb = NULL
                dw = NULL
                
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")
                  p = c(p, p+1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                }
                if("dweibull" %in% lois){
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                }
                
                if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                
                MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])]), 
                               dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])]),
                               dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])]), 
                               ddw), nrow = Kmax)
                
                # print(MGM)
                fuv = matrix(nrow = S-1, ncol = Kmax)
                if("pois" %in% lois){
                  nr = which(lois == "pois")
                  fuv[nr,] = MGM[,1]
                }
                if("geom" %in% lois){
                  nr = which(lois == "geom")
                  fuv[nr,] = MGM[,2]
                }
                if("nbinom" %in% lois){
                  nr = which(lois == "nbinom")
                  fuv[nr,] = MGM[,3]
                }
                if("dweibull" %in% lois){    
                  nr = which(lois == "dweibull")
                  fuv[nr,] = MGM[,4]
                }
                
                -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                   # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                   +  sum(vect.Nijk[i+S*(j-1)+plus]*dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])], log = TRUE))
                   + sum(vect.Nbij[i+S*(j-1)+plus]*log(pgeom(1:Kmax-1, prob = par[(S*(S-2)+g[1])], lower.tail = FALSE)))
                   + sum( Nei[i,]* log(1-puv%*%fuv))
                )
              }
              
              
              ### CONSTRAINTS
              # puv > 0
              u0 = matrix(0, nrow = (S*(S-2))+2*(S-1), ncol = (S*(S-2))+2*(S-1))
              diag(u0) <- 1
              c0 = rep(0,(S*(S-2))+2*(S-1))
              
              # puv < 1
              u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2))+2*(S-1))
              diag(u1) <- -1
              c1 = c(rep(-1,S*(S-2)))
              
              u2 = rbind(u0,u1)
              c2 = c(c0,c1)
              
              ## verification of constraints
              # k2 <- length(c2)
              # n <- dim(u2)[2]
              # for(i in seq_len(k2)) {
              #   j <- which( u2[i,] != 0 )
              #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
              #   cat(" >= " )
              #   cat( c2[i], "\n" )
              # }
              p = c(0,0)
              g = c(0,0)
              nb = c(0,0)
              dw = c(0,0)
              
              init = c()
              if("pois" %in% lois){
                p = 2*which(lois == "pois")-1
                p = c(p, p+1)
                init[p] = c(2, 1)
              }
              if("geom" %in% lois){
                g = 2*which(lois == "geom")-1
                g = c(g, g+1)
                init[g] = c(0.5, 0.4)
              }
              if("nbinom" %in% lois){
                nb = 2*which(lois == "nbinom")-1
                nb = c(nb, nb+1)
                init[nb] = c(2, 4)
              }
              if("dweibull" %in% lois){ 
                dw = 2*which(lois == "dweibull")-1
                dw = c(dw, dw+1)
                init[dw] = c(0.6, 0.8)
              }
              
              ### Log-likelihood optimisation
              CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), init),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
              ### Get the parameters
              parC = CO2$par
              parL = parC[1:(S*(S-2))]
              theta = parC[(S*(S-2)+g[1])]
              
              ### Transformation of parameters
              ## Add parameters of type: pat = 1 - pac - pag
              ## in par vector of optim
              b = 1
              parP = c()
              while( b <= length(parL) ){
                parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
                b = b +S-2
                
              }
              
             ## Add zeros in diagonale 
              c = 1 
              parP0 = c()
              while( c < length(parP) ){
                parP0 = c(parP0, 0,parP[c:(c+(S-1))])
                c = c + S
              }
              parP0 = c(parP0, 0)
              P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
              P[which(is.na(P))] = 0
              # param = array(c(matrix(theta0, nrow = S, ncol = S, byrow = TRUE),matrix(0, nrow = S, ncol = S)), c(S,S,2))
              param[i,j,1] = theta
            }  
            ## No censoring  
          }else{        
            
            ##############
            ## Estimation of the transition matrix
            ##############
            P = Nij/Ni
            P[which(is.na(P))] = 0
            
            ##############
            ## Estimation of theta: 
            #############
            ## vecteur plus for getting all the Nij for each k
            plus = seq(from = 0, to = Kmax*S*S-1, by = S*S )
            ## vect.Nijk is the vector of the array Nijk
            vect.Nijk = as.vector(Nijk)
            ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
            logvrais_theta = function(par){
              -(sum(vect.Nijk[i+S*(j-1)+plus]*dgeom(0:(Kmax-1), prob = par, log = TRUE)))
            }
            
            if (i != j ){ 
              ## Warning for warnings due to par initialization
              # out<-optim(par = 1, logvrais_theta, method = "BFGS")
              x = rep(1:(Kmax), vect.Nijk[i+S*(j-1)+plus])
              param[i,j,1] = 1/(sum(x)/sum(vect.Nijk[i+S*(j-1)+plus]))
              # param[i,j,1] = out$par
            }
            
          } 
          
          
        }else if ( distr[i,j] == "nbinom"){
          
          ## censoring at the end
          if( cens.end == 1 && cens.beg == 0){
            if (S == 2){
              #############################
              ## Estimation of parameters 
              #############################
              plus = seq(from = 0, to = Kmax*S*S-1, by = S*S)
              ## vect.Nijk is the vector of the array Nijk
              vect.Nijk = as.vector(Nijk)
              ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
              ## estimation function of  p_{jt+1v} and theta_{jt+1v}
              Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
              El.diag = seq(from= 1, to = S*S, by = S+1)
              Nij.vect = Nij.vect0[-El.diag]
              # Nij0 = Nij0[c(etat.cens),]
              
              
              ### Separation of parameters for log-likelihood estimation
              Nij1 = Nij.vect[1]
              Nij2 = Nij.vect[2]
              # a = 1
              # while (a < length(Nij.vect)){
              #   Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
              #   Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
              #   a = a + (S-1)
              # }
              # 
              lois = distr[i,]
              lois = lois[-i]
              
             ## Functions to estimate
              logvrais = function(par){
                # m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                # p = cbind(m, 1-rowSums(m))
                # puv = p[i,]
                # # ind1 = seq(from = (S*(S-2))+1, to =(S*(S-2)+5), by=2)
                # ind2 = seq(from = (S*(S-2))+2, to =(S*(S-2)+6), by=2)
                p = NULL
                g = NULL
                nb = NULL
                dw = NULL
                
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                }
                if("dweibull" %in% lois){
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                }
                
                if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                
                MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[p[1]]), 
                               dgeom(0:(Kmax-1), prob = par[g[1]]),
                               dnbinom(0:(Kmax-1), size = par[nb[1]], mu = par[nb[2]]), 
                               ddw), nrow = Kmax)
                
                
                # print(MGM)
                fuv = matrix(nrow = S-1, ncol = Kmax)
                if("pois" %in% lois){
                  nr = which(lois == "pois")
                  fuv[nr,] = MGM[,1]
                }
                if("geom" %in% lois){
                  nr = which(lois == "geom")
                  fuv[nr,] = MGM[,2]
                }
                if("nbinom" %in% lois){
                  nr = which(lois == "nbinom")
                  fuv[nr,] = MGM[,3]
                }
                if("dweibull" %in% lois){    
                  nr = which(lois == "dweibull")
                  fuv[nr,] = MGM[,4]
                }
                
                puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                puv.vect0 = as.vector(t(puv))
                puv.vect = puv.vect0[-El.diag]
                
                -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                   # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                   +  sum(vect.Nijk[i+S*(j-1)+plus]*dnbinom(0:(Kmax-1), size = par[nb[1]], mu = par[nb[2]], log = TRUE))
                   + sum( Nei[i,]* log(1-puv.vect%*%fuv))
                )  
                
                
              }
              
              
              ### CONSTRAINTS
              # # puv > 0
              # u0 = matrix(0, nrow = (S*(S-2))+2*(S-1), ncol = (S*(S-2))+2*(S-1))
              # diag(u0) <- 1
              # c0 = rep(0,(S*(S-2))+2*(S-1))
              # 
              # # puv < 1
              # u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2))+2*(S-1))
              # diag(u1) <- -1
              # c1 = c(rep(-1,S*(S-2)))
              # 
              # u2 = rbind(u0,u1)
              # c2 = c(c0,c1)
              
              ## verification of constraints
              # k2 <- length(c2)
              # n <- dim(u2)[2]
              # for(i in seq_len(k2)) {
              #   j <- which( u2[i,] != 0 )
              #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
              #   cat(" >= " )
              #   cat( c2[i], "\n" )
              # }
              
              
              ### Create the vector for par initialization in optim
              init = c()
              if("pois" %in% lois){
                p = 2*which(lois == "pois")-1
                p = c(p, p+1)
                init[p] = c(2, 1)
              }
              if("geom" %in% lois){
                g = 2*which(lois == "geom")-1
                g = c(g, g+1)
                init[g] = c(0.5, 0.4)
              }
              if("nbinom" %in% lois){
                nb = 2*which(lois == "nbinom")-1
                nb = c(nb, nb+1)
                init[nb] = c(2, 4)
              }
              if("dweibull" %in% lois){ 
                dw = 2*which(lois == "dweibull")-1
                dw = c(dw, dw+1)
                init[dw] = c(0.6, 0.8)
              }
              
              
              ### Log-likelihood optimisation
              CO2<-optim(par = init,fn = logvrais, method ="Nelder-Mead")
              ### Get the parameters
              parC = CO2$par
              alpha = parC[nb[1]]
              theta = parC[nb[2]]
              
              puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
              P = puv
              
              param[i,j,1] = alpha
              param[i,j,2] = theta
              
            }else{
            
              #############################
              ## Estimation of parameters 
              #############################
              plus = seq(from = 0, to = Kmax*S*S-1, by = S*S)
              ## vect.Nijk is the vector of the array Nijk
              vect.Nijk = as.vector(Nijk)
              ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
              ## estimation function of  p_{jt+1v} and theta_{jt+1v}
              Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
              El.diag = seq(from= 1, to = S*S, by = S+1)
              Nij.vect = Nij.vect0[-El.diag]
              # Nij0 = Nij0[c(etat.cens),]
              
              # 
              ### Separation of parameters in order to compute log-likelihood
              Nij1 = c()
              Nij2 = c()
              a = 1
              while (a < length(Nij.vect)){
                Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
                Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
                a = a + (S-1)
              }
              
              lois = distr[i,]
              lois = lois[-i]
              
             ## Functions to estimate
              logvrais = function(par){
                m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                p = cbind(m, 1-rowSums(m))
                puv = p[i,]
                
                p = NULL
                g = NULL
                nb = NULL
                dw = NULL
                
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                }
                if("dweibull" %in% lois){
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                }
                
                if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                
                MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])]), 
                               dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])]),
                               dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])]), 
                               ddw), nrow = Kmax)
                # print(MGM)
                fuv = matrix(nrow = S-1, ncol = Kmax)
                if("pois" %in% lois){
                  nr = which(lois == "pois")
                  fuv[nr,] = MGM[,1]
                }
                if("geom" %in% lois){
                  nr = which(lois == "geom")
                  fuv[nr,] = MGM[,2]
                }
                if("nbinom" %in% lois){
                  nr = which(lois == "nbinom")
                  fuv[nr,] = MGM[,3]
                }
                if("dweibull" %in% lois){    
                  nr = which(lois == "dweibull")
                  fuv[nr,] = MGM[,4]
                }
                -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                   # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                   +  sum(vect.Nijk[i+S*(j-1)+plus]*dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])], log = TRUE))
                   + sum( Nei[i,]* log(1-puv%*%fuv))
                )
              }
              
              
              ### CONSTRAINTS
              # puv > 0
              u0 = matrix(0, nrow = (S*(S-2))+2*(S-1), ncol = (S*(S-2))+2*(S-1))
              diag(u0) <- 1
              c0 = rep(0,(S*(S-2))+2*(S-1))
              
              # puv < 1
              u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2))+2*(S-1))
              diag(u1) <- -1
              c1 = c(rep(-1,S*(S-2)))
              
              u2 = rbind(u0,u1)
              c2 = c(c0,c1)
              
              ## verification of constraints
              # k2 <- length(c2)
              # n <- dim(u2)[2]
              # for(i in seq_len(k2)) {
              #   j <- which( u2[i,] != 0 )
              #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
              #   cat(" >= " )
              #   cat( c2[i], "\n" )
              # }
              p = c(0,0)
              g = c(0,0)
              nb = c(0,0)
              dw = c(0,0)
              
              init = c()
              if("pois" %in% lois){
                p = 2*which(lois == "pois")-1
                p = c(p, p+1)
                init[p] = c(2, 1)
              }
              if("geom" %in% lois){
                g = 2*which(lois == "geom")-1
                g = c(g, g+1)
                init[g] = c(0.5, 0.4)
              }
              if("nbinom" %in% lois){
                nb = 2*which(lois == "nbinom")-1
                nb = c(nb, nb+1)
                init[nb] = c(2, 4)
              }
              if("dweibull" %in% lois){ 
                dw = 2*which(lois == "dweibull")-1
                dw = c(dw, dw+1)
                init[dw] = c(0.6, 0.8)
              }
              
              ### Log-likelihood optimisation
              CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), init),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
              ### Get the parameters
              parC = CO2$par
              parL = parC[1:(S*(S-2))]
              alpha = parC[(S*(S-2)+nb[1])]
              theta = parC[(S*(S-2)+nb[2])]
              
              ### Transformation of parameters
                ## Add parameters of type: pat = 1 - pac - pag
                ## in par vector of optim
              b = 1
              parP = c()
              while( b <= length(parL) ){
                parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
                b = b +S-2
                
              }
              
             ## Add zeros in diagonale 
              c = 1 
              parP0 = c()
              while( c < length(parP) ){
                parP0 = c(parP0, 0,parP[c:(c+(S-1))])
                c = c + S
              }
              parP0 = c(parP0, 0)
              P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
              P[which(is.na(P))] = 0
              # param = array(c(matrix(theta0, nrow = S, ncol = S, byrow = TRUE),matrix(0, nrow = S, ncol = S)), c(S,S,2))
              param[i,j,1] = alpha
              param[i,j,2] = theta
              
            }
            ## censoring at the beginning 
          }else if( cens.end == 0 && cens.beg == 1 ){   
              
              ########################
              ## Estimation of the transition matrix
              #######################
              P = Nij/Ni
              P[which(is.na(P))] = 0
              
              plus = seq(from = 0, to = Kmax*S*S-1, by = S*S )
              ## vect.Nijk is the vector of the array Nijk
              vect.Nijk = as.vector(Nijk)
              ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
              ## In the same way, we create the vector  vect.Nbij : 
              ## vect.Nbij is the vector of the array Nbij
              vect.Nbij = as.vector((Nbij))
              ## the Nbij[i,j,] correspond for the vector Nbij to vect.Nbij[i+S*(j-1)+plus]
              ########################
              ## Estimation of parameters  of the distribution 
              ########################
              logvrais_theta = function(par){
                -(sum(vect.Nijk[i+S*(j-1)+plus]*dnbinom(0:(Kmax-1), size = par[1], mu = par[2], log = TRUE))
                  + sum(vect.Nbij[i+S*(j-1)+plus]*log(pnbinom(1:Kmax-1, size = par[1], mu = par[2], lower.tail = FALSE)))
                )
              }
              
              out<-optim(par = c(5,2), logvrais_theta)
              param[i,j,1] = out$par[1]
              param[i,j,2] = out$par[2]
              
              ## censoring at the end and at the beginning
            }else if( cens.end == 1 && cens.beg == 1 ){
              
              if (S == 2){
                #############################
                ## Estimation of parameters 
                #############################
                plus = seq(from = 0, to = Kmax*S*S-1, by = S*S)
                ## vect.Nijk is the vector of the array Nijk
                vect.Nijk = as.vector(Nijk)
                ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
                ## estimation function of  p_{jt+1v} and theta_{jt+1v}
                Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
                El.diag = seq(from= 1, to = S*S, by = S+1)
                Nij.vect = Nij.vect0[-El.diag]
                # Nij0 = Nij0[c(etat.cens),]
                ## vect.Nbij is the vector of the array Nbij
                vect.Nbij = as.vector((Nbij))
                ## the Nbij[i,j,] correspond for the vector Nbij to vect.Nbij[i+S*(j-1)+plus]
                
                ### Separation of parameters in order to compute the log-likelihood
                Nij1 = Nij.vect[1]
                Nij2 = Nij.vect[2]
                # a = 1
                # while (a < length(Nij.vect)){
                #   Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
                #   Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
                #   a = a + (S-1)
                # }
                # 
                lois = distr[i,]
                lois = lois[-i]
                
               ## Functions to estimate
                logvrais = function(par){
                  # m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                  # p = cbind(m, 1-rowSums(m))
                  # puv = p[i,]
                  # # ind1 = seq(from = (S*(S-2))+1, to =(S*(S-2)+5), by=2)
                  # ind2 = seq(from = (S*(S-2))+2, to =(S*(S-2)+6), by=2)
                  p = NULL
                  g = NULL
                  nb = NULL
                  dw = NULL
                  
                  if("pois" %in% lois){
                    p = 2*which(lois == "pois")-1
                    p = c(p, p+1)
                  }
                  if("geom" %in% lois){
                    g = 2*which(lois == "geom")-1
                    g = c(g, g+1)
                  }
                  if("nbinom" %in% lois){
                    nb = 2*which(lois == "nbinom")-1
                    nb = c(nb, nb+1)
                  }
                  if("dweibull" %in% lois){
                    dw = 2*which(lois == "dweibull")-1
                    dw = c(dw, dw+1)
                  }
                  
                  if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                  
                  MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[p[1]]), 
                                 dgeom(0:(Kmax-1), prob = par[g[1]]),
                                 dnbinom(0:(Kmax-1), size = par[nb[1]], mu = par[nb[2]]), 
                                 ddw), nrow = Kmax)
                  
                  
                  # print(MGM)
                  fuv = matrix(nrow = S-1, ncol = Kmax)
                  if("pois" %in% lois){
                    nr = which(lois == "pois")
                    fuv[nr,] = MGM[,1]
                  }
                  if("geom" %in% lois){
                    nr = which(lois == "geom")
                    fuv[nr,] = MGM[,2]
                  }
                  if("nbinom" %in% lois){
                    nr = which(lois == "nbinom")
                    fuv[nr,] = MGM[,3]
                  }
                  if("dweibull" %in% lois){    
                    nr = which(lois == "dweibull")
                    fuv[nr,] = MGM[,4]
                  }
                  
                  puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                  puv.vect0 = as.vector(t(puv))
                  puv.vect = puv.vect0[-El.diag]
                  
                  -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                     # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                     +  sum(vect.Nijk[i+S*(j-1)+plus]*dnbinom(0:(Kmax-1), size = par[nb[1]], mu = par[nb[2]], log = TRUE))
                     + sum(vect.Nbij[i+S*(j-1)+plus]*log(pnbinom(1:Kmax-1, size = par[nb[1]], mu = par[nb[2]], lower.tail = FALSE)))
                     + sum( Nei[i,]* log(1-puv.vect%*%fuv))
                  )  
                  
                  
                }
                
                
                ### CONSTRAINTS
                # # puv > 0
                # u0 = matrix(0, nrow = (S*(S-2))+2*(S-1), ncol = (S*(S-2))+2*(S-1))
                # diag(u0) <- 1
                # c0 = rep(0,(S*(S-2))+2*(S-1))
                # 
                # # puv < 1
                # u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2))+2*(S-1))
                # diag(u1) <- -1
                # c1 = c(rep(-1,S*(S-2)))
                # 
                # u2 = rbind(u0,u1)
                # c2 = c(c0,c1)
                
                ## verification of constraints
                # k2 <- length(c2)
                # n <- dim(u2)[2]
                # for(i in seq_len(k2)) {
                #   j <- which( u2[i,] != 0 )
                #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
                #   cat(" >= " )
                #   cat( c2[i], "\n" )
                # }
                
                
                ### Vector to initialize par in optim
                init = c()
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                  init[p] = c(2, 1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                  init[g] = c(0.5, 0.4)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                  init[nb] = c(2, 4)
                }
                if("dweibull" %in% lois){ 
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                  init[dw] = c(0.6, 0.8)
                }
                
                
                ### Log-likelihood optimisation
                CO2<-optim(par = init,fn = logvrais, method ="Nelder-Mead")
                ### Get the parameters
                parC = CO2$par
                alpha = parC[nb[1]]
                theta = parC[nb[2]]
                
                puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                P = puv
                
                param[i,j,1] = alpha
                param[i,j,2] = theta
                
                
                
              }else{
                #############################
                ## Estimation of parameters 
                #############################
                plus = seq(from = 0, to = Kmax*S*S-1, by = S*S )
                ## vect.Nijk is the vector of the array Nijk
                vect.Nijk = as.vector(Nijk)
                ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
                ## In the same way, we create the vector  vect.Nbij : 
                ## vect.Nbij is the vector of the array Nbij
                vect.Nbij = as.vector((Nbij))
                ## the Nbij[i,j,] correspond for the vector Nbij to vect.Nbij[i+S*(j-1)+plus]
                Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
                El.diag = seq(from= 1, to = S*S, by = S+1)
                Nij.vect = Nij.vect0[-El.diag]
                # Nij0 = Nij0[c(etat.cens),]
                
                # 
                ### Separation of parameters in order to compute log-likelihood
                Nij1 = c()
                Nij2 = c()
                a = 1
                while (a < length(Nij.vect)){
                  Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
                  Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
                  a = a + (S-1)
                }
                
                lois = distr[i,]
                lois = lois[-i]
                
               ## Functions to estimate
                logvrais = function(par){
                  m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                  p = cbind(m, 1-rowSums(m))
                  puv = p[i,]
                  p = NULL
                  g = NULL
                  nb = NULL
                  dw = NULL
                  
                  if("pois" %in% lois){
                    p = 2*which(lois == "pois")-1
                    p = c(p, p+1)
                  }
                  if("geom" %in% lois){
                    g = 2*which(lois == "geom")-1
                    g = c(g, g+1)
                  }
                  if("nbinom" %in% lois){
                    nb = 2*which(lois == "nbinom")-1
                    nb = c(nb, nb+1)
                  }
                  if("dweibull" %in% lois){
                    dw = 2*which(lois == "dweibull")-1
                    dw = c(dw, dw+1)
                  }
                  
                  if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                  
                  MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])]), 
                                 dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])]),
                                 dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])]), 
                                 ddw), nrow = Kmax)
                  
                  # print(MGM)
                  fuv = matrix(nrow = S-1, ncol = Kmax)
                  if("pois" %in% lois){
                    nr = which(lois == "pois")
                    fuv[nr,] = MGM[,1]
                  }
                  if("geom" %in% lois){
                    nr = which(lois == "geom")
                    fuv[nr,] = MGM[,2]
                  }
                  if("nbinom" %in% lois){
                    nr = which(lois == "nbinom")
                    fuv[nr,] = MGM[,3]
                  }
                  if("dweibull" %in% lois){    
                    nr = which(lois == "dweibull")
                    fuv[nr,] = MGM[,4]
                  }
                  -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                     # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                     + sum(vect.Nijk[i+S*(j-1)+plus]*dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])], log = TRUE))
                     + sum(vect.Nbij[i+S*(j-1)+plus]*log(pnbinom(1:Kmax-1, size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])], lower.tail = FALSE)))
                     + sum( Nei[i,]* log(1-puv%*%fuv))
                  )
                }
                
                
                ### CONSTRAINTS
                # puv > 0
                u0 = matrix(0, nrow = (S*(S-2))+2*(S-1), ncol = (S*(S-2))+2*(S-1))
                diag(u0) <- 1
                c0 = rep(0,(S*(S-2))+2*(S-1))
                
                # puv < 1
                u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2))+2*(S-1))
                diag(u1) <- -1
                c1 = c(rep(-1,S*(S-2)))
                
                u2 = rbind(u0,u1)
                c2 = c(c0,c1)
                
                ## verification of constraints
                # k2 <- length(c2)
                # n <- dim(u2)[2]
                # for(i in seq_len(k2)) {
                #   j <- which( u2[i,] != 0 )
                #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
                #   cat(" >= " )
                #   cat( c2[i], "\n" )
                # }
                
                p = c(0,0)
                g = c(0,0)
                nb = c(0,0)
                dw = c(0,0)
                init = c()
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                  init[p] = c(2, 1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                  init[g] = c(0.5, 0.4)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                  init[nb] = c(2, 4)
                }
                if("dweibull" %in% lois){ 
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                  init[dw] = c(0.6, 0.8)
                }
                
                ### Log-likelihood optimisation
                CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), init),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
                ### Get the parameters
                parC = CO2$par
                parL = parC[1:(S*(S-2))]
                alpha = parC[(S*(S-2)+nb[1])]
                theta = parC[(S*(S-2)+nb[2])]
                
                ### Transformation of parameters
                ## Add parameters of type: pat = 1 - pac - pag
                ## in par vector of optim
                b = 1
                parP = c()
                while( b <= length(parL) ){
                  parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
                  b = b +S-2
                  
                }
                
               ## Add zeros on diagonal
                c = 1 
                parP0 = c()
                while( c < length(parP) ){
                  parP0 = c(parP0, 0,parP[c:(c+(S-1))])
                  c = c + S
                }
                parP0 = c(parP0, 0)
                P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
                P[which(is.na(P))] = 0
                # param = array(c(matrix(theta0, nrow = S, ncol = S, byrow = TRUE),matrix(0, nrow = S, ncol = S)), c(S,S,2))
                param[i,j,1] = alpha
                param[i,j,2] = theta
                
              }
              
              ## No censoring  
            }else{        
              
              ##############
              ## Estimation of the transition matrix
              ##############
              P = Nij/Ni
              P[which(is.na(P))] = 0
              
              ##############
              ## Estimation of theta
              #############
              ## vecteur plus for getting all the Nij for each k
              plus = seq(from = 0, to = Kmax*S*S-1, by = S*S )
              ## vect.Nijk is the vector of the array Nijk
              vect.Nijk = as.vector(Nijk)
              ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
              logvrais_theta = function(par){
                -(sum(vect.Nijk[i+S*(j-1)+plus]*dnbinom(0:(Kmax-1), size = par[1], mu = par[2], log = TRUE)))
              }
              
              if (i != j ){ 
                # Warning for warnings due to par initialization
                out<-optim(par = c(2,5), logvrais_theta)
                param[i,j,1] = out$par[1]
                param[i,j,2] = out$par[2]
              }
              
            } 
            
            
          }else if ( distr[i,j] == "dweibull" ){
            
            ## censoring at the end
            if( cens.end == 1 && cens.beg == 0){
              if( S == 2 ){
                #############################
                ## Estimation of parameters 
                #############################
                plus = seq(from = 0, to = Kmax*S*S-1, by = S*S)
                ## vect.Nijk is the vector of the array Nijk
                vect.Nijk = as.vector(Nijk)
                ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
                ## estimation function of  p_{jt+1v} and theta_{jt+1v}
                Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
                El.diag = seq(from= 1, to = S*S, by = S+1)
                Nij.vect = Nij.vect0[-El.diag]
                # Nij0 = Nij0[c(etat.cens),]
                
                
                ### Separation of parameters in order to compute log-likelihood
                Nij1 = Nij.vect[1]
                Nij2 = Nij.vect[2]
                # a = 1
                # while (a < length(Nij.vect)){
                #   Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
                #   Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
                #   a = a + (S-1)
                # }
                # 
                lois = distr[i,]
                lois = lois[-i]
                
               ## Functions to estimate
                logvrais = function(par){
                  # m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                  # p = cbind(m, 1-rowSums(m))
                  # puv = p[i,]
                  # # ind1 = seq(from = (S*(S-2))+1, to =(S*(S-2)+5), by=2)
                  # ind2 = seq(from = (S*(S-2))+2, to =(S*(S-2)+6), by=2)
                  p = NULL
                  g = NULL
                  nb = NULL
                  dw = NULL
                  
                  if("pois" %in% lois){
                    p = 2*which(lois == "pois")-1
                    p = c(p, p+1)
                  }
                  if("geom" %in% lois){
                    g = 2*which(lois == "geom")-1
                    g = c(g, g+1)
                  }
                  if("nbinom" %in% lois){
                    nb = 2*which(lois == "nbinom")-1
                    nb = c(nb, nb+1)
                  }
                  if("dweibull" %in% lois){
                    dw = 2*which(lois == "dweibull")-1
                    dw = c(dw, dw+1)
                  }
                  
                  if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                  
                  MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[p[1]]), 
                                 dgeom(0:(Kmax-1), prob = par[g[1]]),
                                 dnbinom(0:(Kmax-1), size = par[nb[1]], mu = par[nb[2]]), 
                                 ddw), nrow = Kmax)
                  
                  
                  # print(MGM)
                  fuv = matrix(nrow = S-1, ncol = Kmax)
                  if("pois" %in% lois){
                    nr = which(lois == "pois")
                    fuv[nr,] = MGM[,1]
                  }
                  if("geom" %in% lois){
                    nr = which(lois == "geom")
                    fuv[nr,] = MGM[,2]
                  }
                  if("nbinom" %in% lois){
                    nr = which(lois == "nbinom")
                    fuv[nr,] = MGM[,3]
                  }
                  if("dweibull" %in% lois){    
                    nr = which(lois == "dweibull")
                    fuv[nr,] = MGM[,4]
                  }
                  
                  puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                  puv.vect0 = as.vector(t(puv))
                  puv.vect = puv.vect0[-El.diag]
                  
                  -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                     # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                     +  sum(vect.Nijk[i+S*(j-1)+plus]*log(ddweibull(0:(Kmax-1), q = par[dw[1]], beta = par[dw[2]])))
                     + sum( Nei[i,]* log(1-puv.vect%*%fuv))
                  )  
                  
                  
                }
                
                
                ### CONSTRAINTS
                # # puv > 0
                # u0 = matrix(0, nrow = (S*(S-2))+2*(S-1), ncol = (S*(S-2))+2*(S-1))
                # diag(u0) <- 1
                # c0 = rep(0,(S*(S-2))+2*(S-1))
                # 
                # # puv < 1
                # u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2))+2*(S-1))
                # diag(u1) <- -1
                # c1 = c(rep(-1,S*(S-2)))
                # 
                # u2 = rbind(u0,u1)
                # c2 = c(c0,c1)
                
                ## verification of constraints
                # k2 <- length(c2)
                # n <- dim(u2)[2]
                # for(i in seq_len(k2)) {
                #   j <- which( u2[i,] != 0 )
                #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
                #   cat(" >= " )
                #   cat( c2[i], "\n" )
                # }
                
                
                ### Vector to initialize par in optim
                init = c()
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                  init[p] = c(2, 1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                  init[g] = c(0.5, 0.4)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                  init[nb] = c(2, 4)
                }
                if("dweibull" %in% lois){ 
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                  init[dw] = c(0.6, 0.8)
                }
                
                
                ### Log-likelihood optimisation
                CO2<-optim(par = init,fn = logvrais, method ="Nelder-Mead")
                ### Get the parameters
                parC = CO2$par
                q = parC[dw[1]]
                beta = parC[dw[2]]
                
                param[i,j,1] = q
                param[i,j,2] = beta
                puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                P = puv
                
              }else{
                #############################
                ## Estimation of parameters 
                #############################
                # tNijk = aperm(Nijk, c(2,1,3))
                # Nijk.vect = as.vector(tNijk)
                # n = 1
                # d = 1
                # Diag = c()
                # while ( n < S*S*Kmax ){
                #   Diag = c(Diag, seq(from = n, to = d*S*S, by = S+1))
                #   d=d+1
                #   n=n+S*S
                # }
                # Nijk.vect = Nijk.vect[-Diag]
                plus = seq(from = 0, to = Kmax*S*S-1, by = S*S)
                ## vect.Nijk is the vector of the array Nijk
                vect.Nijk = as.vector(Nijk)
                ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
                ## estimation function of  p_{jt+1v} and theta_{jt+1v}
                Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
                El.diag = seq(from= 1, to = S*S, by = S+1)
                Nij.vect = Nij.vect0[-El.diag]
                # Nij0 = Nij0[c(etat.cens),]
                
                # 
                ### Separation of parameters in order to compute the log-likelihood
                Nij1 = c()
                Nij2 = c()
                a = 1
                while (a < length(Nij.vect)){
                  Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
                  Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
                  a = a + (S-1)
                }
                
                lois = distr[i,]
                lois = lois[-i]
                
               ## Functions to estimate
                logvrais = function(par){
                  m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                  p = cbind(m, 1-rowSums(m))
                  puv = p[i,]
                  p = NULL
                  g = NULL
                  nb = NULL
                  dw = NULL
                  
                  if("pois" %in% lois){
                    p = 2*which(lois == "pois")-1
                    p = c(p, p+1)
                  }
                  if("geom" %in% lois){
                    g = 2*which(lois == "geom")-1
                    g = c(g, g+1)
                  }
                  if("nbinom" %in% lois){
                    nb = 2*which(lois == "nbinom")-1
                    nb = c(nb, nb+1)
                  }
                  if("dweibull" %in% lois){
                    dw = 2*which(lois == "dweibull")-1
                    dw = c(dw, dw+1)
                  }
                  
                  if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                  
                  MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])]), 
                                 dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])]),
                                 dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])]), 
                                 ddw), nrow = Kmax)
                  
                  # print(MGM)
                  fuv = matrix(nrow = S-1, ncol = Kmax)
                  if("pois" %in% lois){
                    nr = which(lois == "pois")
                    fuv[nr,] = MGM[,1]
                  }
                  if("geom" %in% lois){
                    nr = which(lois == "geom")
                    fuv[nr,] = MGM[,2]
                  }
                  if("nbinom" %in% lois){
                    nr = which(lois == "nbinom")
                    fuv[nr,] = MGM[,3]
                  }
                  if("dweibull" %in% lois){    
                    nr = which(lois == "dweibull")
                    fuv[nr,] = MGM[,4]
                  }
                  -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                     # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                     +  sum(vect.Nijk[i+S*(j-1)+plus]*log(ddweibull(1:(Kmax), q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])))
                     + sum( Nei[i,]* log(1-puv%*%fuv))
                  )
                }
                
                
                ### CONSTRAINTS
                # puv > 0
                u0 = matrix(0, nrow = (S*(S-2))+2*(S-1), ncol = (S*(S-2))+2*(S-1))
                diag(u0) <- 1
                c0 = rep(0,(S*(S-2))+2*(S-1))
                
                # puv < 1
                u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2))+2*(S-1))
                diag(u1) <- -1
                c1 = c(rep(-1,S*(S-2)))
                
                u2 = rbind(u0,u1)
                c2 = c(c0,c1)
                
                ## verification of constraints
                # k2 <- length(c2)
                # n <- dim(u2)[2]
                # for(i in seq_len(k2)) {
                #   j <- which( u2[i,] != 0 )
                #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
                #   cat(" >= " )
                #   cat( c2[i], "\n" )
                # }
                
                p = c(0,0)
                g = c(0,0)
                nb = c(0,0)
                dw = c(0,0)
                init = c()
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                  init[p] = c(2, 1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                  init[g] = c(0.5, 0.4)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                  init[nb] = c(2, 4)
                }
                if("dweibull" %in% lois){ 
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                  init[dw] = c(0.6, 0.8)
                }
                
                ### Log-likelihood optimisation
                CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), init),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
                ### Get the parameters
                parC = CO2$par
                parL = parC[1:(S*(S-2))]
                q = parC[(S*(S-2)+dw[1])]
                beta = parC[(S*(S-2)+dw[2])]
                
                ### Transformation of parameters
                ## Add parameters of type: pat = 1 - pac - pag
                ## in par vector of optim
                b = 1
                parP = c()
                while( b <= length(parL) ){
                  parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
                  b = b +S-2
                  
                }
                
               ## Add zeros in diagonale 
                c = 1 
                parP0 = c()
                while( c < length(parP) ){
                  parP0 = c(parP0, 0,parP[c:(c+(S-1))])
                  c = c + S
                }
                parP0 = c(parP0, 0)
                P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
                P[which(is.na(P))] = 0
                # param = array(c(matrix(theta0, nrow = S, ncol = S, byrow = TRUE),matrix(0, nrow = S, ncol = S)), c(S,S,2))
                param[i,j,1] = q
                param[i,j,2] = beta
                
              }
              
              ## censoring at the beginning 
            }else if( cens.end == 0 && cens.beg == 1 ){   
              
              ########################
              ## Estimation of the transition matrix
              #######################
              P = Nij/Ni
              P[which(is.na(P))] = 0
              
              plus = seq(from = 0, to = Kmax*S*S-1, by = S*S )
              ## vect.Nijk is the vector of the array Nijk
              vect.Nijk = as.vector(Nijk)
              ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
              ## In the same way, we create the vector  vect.Nbij : 
              ## vect.Nbij is the vector of the array Nbij
              vect.Nbij = as.vector((Nbij))
              ## the Nbij[i,j,] correspond for the vector Nbij to vect.Nbij[i+S*(j-1)+plus]
              ########################
              ## Estimation of parameters  of the distribution 
              ########################
              logvrais_theta = function(par){
                -(sum(vect.Nijk[i+S*(j-1)+plus]*log(ddweibull(1:(Kmax), q = par[1], beta = par[2])))
                  + sum(vect.Nbij[i+S*(j-1)+plus]*log(1-pdweibull(1:Kmax, q = par[1], beta = par[2])))
                )
              }
              
              out<-optim(par = c(0.3,0.6), logvrais_theta)
              param[i,j,1] = out$par[1]
              param[i,j,2] = out$par[2]
              
              ## censoring at the end and at the beginning
            }else if( cens.end == 1 && cens.beg == 1 ){
              
              if (S == 2){
                #############################
                ## Estimation of parameters 
                #############################
                plus = seq(from = 0, to = Kmax*S*S-1, by = S*S)
                ## vect.Nijk is the vector of the array Nijk
                vect.Nijk = as.vector(Nijk)
                ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
                ## estimation function of  p_{jt+1v} and theta_{jt+1v}
                Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
                El.diag = seq(from= 1, to = S*S, by = S+1)
                Nij.vect = Nij.vect0[-El.diag]
                # Nij0 = Nij0[c(etat.cens),]
                ## vect.Nbij is the vector of the array Nbij
                vect.Nbij = as.vector((Nbij))
                ## the Nbij[i,j,] correspond for the vector Nbij to vect.Nbij[i+S*(j-1)+plus]
                
                
                ### Separation of parameters in order to compute log-likelihood
                Nij1 = Nij.vect[1]
                Nij2 = Nij.vect[2]
                # a = 1
                # while (a < length(Nij.vect)){
                #   Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
                #   Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
                #   a = a + (S-1)
                # }
                # 
                lois = distr[i,]
                lois = lois[-i]
                
               ## Functions to estimate
                logvrais = function(par){
                  # m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                  # p = cbind(m, 1-rowSums(m))
                  # puv = p[i,]
                  # # ind1 = seq(from = (S*(S-2))+1, to =(S*(S-2)+5), by=2)
                  # ind2 = seq(from = (S*(S-2))+2, to =(S*(S-2)+6), by=2)
                  p = NULL
                  g = NULL
                  nb = NULL
                  dw = NULL
                  
                  if("pois" %in% lois){
                    p = 2*which(lois == "pois")-1
                    p = c(p, p+1)
                  }
                  if("geom" %in% lois){
                    g = 2*which(lois == "geom")-1
                    g = c(g, g+1)
                  }
                  if("nbinom" %in% lois){
                    nb = 2*which(lois == "nbinom")-1
                    nb = c(nb, nb+1)
                  }
                  if("dweibull" %in% lois){
                    dw = 2*which(lois == "dweibull")-1
                    dw = c(dw, dw+1)
                  }
                  
                  if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                  
                  MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[p[1]]), 
                                 dgeom(0:(Kmax-1), prob = par[g[1]]),
                                 dnbinom(0:(Kmax-1), size = par[nb[1]], mu = par[nb[2]]), 
                                 ddw), nrow = Kmax)
                  
                  
                  # print(MGM)
                  fuv = matrix(nrow = S-1, ncol = Kmax)
                  if("pois" %in% lois){
                    nr = which(lois == "pois")
                    fuv[nr,] = MGM[,1]
                  }
                  if("geom" %in% lois){
                    nr = which(lois == "geom")
                    fuv[nr,] = MGM[,2]
                  }
                  if("nbinom" %in% lois){
                    nr = which(lois == "nbinom")
                    fuv[nr,] = MGM[,3]
                  }
                  if("dweibull" %in% lois){    
                    nr = which(lois == "dweibull")
                    fuv[nr,] = MGM[,4]
                  }
                  
                  puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                  puv.vect0 = as.vector(t(puv))
                  puv.vect = puv.vect0[-El.diag]
                  
                  -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                     # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                     +  sum(vect.Nijk[i+S*(j-1)+plus]*log(ddweibull(0:(Kmax-1), q = par[dw[1]], beta = par[dw[2]])))
                     + sum(vect.Nbij[i+S*(j-1)+plus]*log(1-pdweibull(1:Kmax, q = par[dw[1]], beta = par[dw[2]])))
                     + sum( Nei[i,]* log(1-puv.vect%*%fuv))
                  )  
                  
                  
                }
                
                
                ### CONSTRAINTS
                # # puv > 0
                # u0 = matrix(0, nrow = (S*(S-2))+2*(S-1), ncol = (S*(S-2))+2*(S-1))
                # diag(u0) <- 1
                # c0 = rep(0,(S*(S-2))+2*(S-1))
                # 
                # # puv < 1
                # u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2))+2*(S-1))
                # diag(u1) <- -1
                # c1 = c(rep(-1,S*(S-2)))
                # 
                # u2 = rbind(u0,u1)
                # c2 = c(c0,c1)
                
                ## verification of constraints
                # k2 <- length(c2)
                # n <- dim(u2)[2]
                # for(i in seq_len(k2)) {
                #   j <- which( u2[i,] != 0 )
                #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
                #   cat(" >= " )
                #   cat( c2[i], "\n" )
                # }
                
                
                ### Vector to initialize par in optim
                init = c()
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                  init[p] = c(2, 1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                  init[g] = c(0.5, 0.4)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                  init[nb] = c(2, 4)
                }
                if("dweibull" %in% lois){ 
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                  init[dw] = c(0.6, 0.8)
                }
                
                
                ### Log-likelihood optimisation
                CO2<-optim(par = init,fn = logvrais, method ="Nelder-Mead")
                ### Get the parameters
                parC = CO2$par
                q = parC[dw[1]]
                beta = parC[dw[2]]
                
                puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
                P = puv
                
                param[i,j,1] = q
                param[i,j,2] = beta
                
              }else{
                #############################
                ## Estimation of parameters 
                #############################
                # tNijk = aperm(Nijk, c(2,1,3))
                # Nijk.vect = as.vector(tNijk)
                # n = 1
                # d = 1
                # Diag = c()
                # while ( n < S*S*Kmax ){
                #   Diag = c(Diag, seq(from = n, to = d*S*S, by = S+1))
                #   d=d+1
                #   n=n+S*S
                # }
                # Nijk.vect = Nijk.vect[-Diag]
                plus = seq(from = 0, to = Kmax*S*S-1, by = S*S )
                ## vect.Nijk is the vector of the array Nijk
                vect.Nijk = as.vector(Nijk)
                ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
                ## In the same way, we create the vector  vect.Nbij : 
                ## vect.Nbij is the vector of the array Nbij
                vect.Nbij = as.vector((Nbij))
                ## the Nbij[i,j,] correspond for the vector Nbij to vect.Nbij[i+S*(j-1)+plus]
                Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
                El.diag = seq(from= 1, to = S*S, by = S+1)
                Nij.vect = Nij.vect0[-El.diag]
                # Nij0 = Nij0[c(etat.cens),]
                
                ### Separation of parameters puv for the likelihood estimation
                Nij1 = c()
                Nij2 = c()
                a = 1
                while (a < length(Nij.vect)){
                  Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
                  Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
                  a = a + (S-1)
                }
                
                lois = distr[i,]
                lois = lois[-i]
                
               ## Functions to estimate
                logvrais = function(par){
                  m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
                  p = cbind(m, 1-rowSums(m))
                  puv = p[i,]
                  p = NULL
                  g = NULL
                  nb = NULL
                  dw = NULL
                  
                  
                  if("pois" %in% lois){
                    p = 2*which(lois == "pois")-1
                    p = c(p, p+1)
                  }
                  if("geom" %in% lois){
                    g = 2*which(lois == "geom")-1
                    g = c(g, g+1)
                  }
                  if("nbinom" %in% lois){
                    nb = 2*which(lois == "nbinom")-1
                    nb = c(nb, nb+1)
                  }
                  if("dweibull" %in% lois){
                    dw = 2*which(lois == "dweibull")-1
                    dw = c(dw, dw+1)
                  }
                  if(is.null(dw)){ddw = rep(0, Kmax)}else{ddw = ddweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])}
                  
                  MGM = matrix(c(dpois(0:(Kmax-1), lambda = par[(S*(S-2)+p[1])]), 
                                 dgeom(0:(Kmax-1), prob = par[(S*(S-2)+g[1])]),
                                 dnbinom(0:(Kmax-1), size = par[(S*(S-2)+nb[1])], mu = par[(S*(S-2)+nb[2])]), 
                                 ddw), nrow = Kmax)
                  # print(MGM)
                  fuv = matrix(nrow = S-1, ncol = Kmax)
                  if("pois" %in% lois){
                    nr = which(lois == "pois")
                    fuv[nr,] = MGM[,1]
                  }
                  if("geom" %in% lois){
                    nr = which(lois == "geom")
                    fuv[nr,] = MGM[,2]
                  }
                  if("nbinom" %in% lois){
                    nr = which(lois == "nbinom")
                    fuv[nr,] = MGM[,3]
                  }
                  if("dweibull" %in% lois){    
                    nr = which(lois == "dweibull")
                    fuv[nr,] = MGM[,4]
                  }
                  -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                     # + sum(Nijk.vect * dpois(rep(0:(Kmax-1), each = (S-1)*S), lambda = par[(S*(S-2)+1)], log = TRUE))
                     + sum(vect.Nijk[i+S*(j-1)+plus]*log(ddweibull(1:(Kmax), q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])))
                     + sum(vect.Nbij[i+S*(j-1)+plus]*log(1-pdweibull(1:Kmax, q = par[(S*(S-2)+dw[1])], beta = par[(S*(S-2)+dw[2])])))
                     + sum( Nei[i,]* log(1-puv%*%fuv))
                  )
                }
                
                
                ### CONSTRAINTS
                # puv > 0
                u0 = matrix(0, nrow = (S*(S-2))+2*(S-1), ncol = (S*(S-2))+2*(S-1))
                diag(u0) <- 1
                c0 = rep(0,(S*(S-2))+2*(S-1))
                
                # puv < 1
                u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2))+2*(S-1))
                diag(u1) <- -1
                c1 = c(rep(-1,S*(S-2)))
                
                u2 = rbind(u0,u1)
                c2 = c(c0,c1)
                
                ## verification of constraints
                # k2 <- length(c2)
                # n <- dim(u2)[2]
                # for(i in seq_len(k2)) {
                #   j <- which( u2[i,] != 0 )
                #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
                #   cat(" >= " )
                #   cat( c2[i], "\n" )
                # }
                p = c(0,0)
                g = c(0,0)
                nb = c(0,0)
                dw = c(0,0)
                
                init = c()
                if("pois" %in% lois){
                  p = 2*which(lois == "pois")-1
                  p = c(p, p+1)
                  init[p] = c(2, 1)
                }
                if("geom" %in% lois){
                  g = 2*which(lois == "geom")-1
                  g = c(g, g+1)
                  init[g] = c(0.5, 0.4)
                }
                if("nbinom" %in% lois){
                  nb = 2*which(lois == "nbinom")-1
                  nb = c(nb, nb+1)
                  init[nb] = c(2, 4)
                }
                if("dweibull" %in% lois){ 
                  dw = 2*which(lois == "dweibull")-1
                  dw = c(dw, dw+1)
                  init[dw] = c(0.6, 0.8)
                }
                
                ### Log-likelihood optimisation
                CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), init),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
                ### Get the parameters
                parC = CO2$par
                parL = parC[1:(S*(S-2))]
                q = parC[(S*(S-2)+dw[1])]
                beta = parC[(S*(S-2)+dw[2])]
                
                ### Transformation of parameters
                ## Add parameters of type: pat = 1 - pac - pag
                ## in par vector of optim
                b = 1
                parP = c()
                while( b <= length(parL) ){
                  parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
                  b = b +S-2
                  
                }
                
               ## Add zeros on diagonal 
                c = 1 
                parP0 = c()
                while( c < length(parP) ){
                  parP0 = c(parP0, 0,parP[c:(c+(S-1))])
                  c = c + S
                }
                parP0 = c(parP0, 0)
                P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
                P[which(is.na(P))] = 0
                # param = array(c(matrix(theta0, nrow = S, ncol = S, byrow = TRUE),matrix(0, nrow = S, ncol = S)), c(S,S,2))
                param[i,j,1] = q
                param[i,j,2] = beta
                
              }
              
              ## No censoring  
            }else{        
              
              ##############
              ## Estimation of the transition matrix
              ##############
              P = Nij/Ni
              P[which(is.na(P))] = 0
              
              ##############
              ## Estimation of theta 
              #############
              ## vecteur plus for getting all the Nij for each k
              plus = seq(from = 0, to = Kmax*S*S-1, by = S*S )
              ## vect.Nijk is the vector of the array Nijk
              vect.Nijk = as.vector(Nijk)
              ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
              logvrais_theta = function(par){
                -(sum(vect.Nijk[i+S*(j-1)+plus]*log(ddweibull(1:(Kmax), q = par[1], beta = par[2]))))
              }
              
              if (i != j ){ 
                # Warning for warnings due to par initialization
                out<-optim(par = c(0.2,0.5), logvrais_theta)
                param[i,j,1] = out$par[1]
                param[i,j,2] = out$par[2]
              }
              
            }
            
            
          }else{
            stop("The name of the distribution is not valid")
          }
          
        }
        
      }
    }
    
    
    ## semi-markovian kernel
    kernel = .kernel_param_fij(distr = distr, param = param, Kmax = Kmax, pij = P, S = S)
    q = kernel$q
    ## Inital law
    if(nbSeq > 1){
      init = Nstart/sum(Nstart)
    }else{
      ## computation of the limit  distribution
      init = .limit.law(q = q, pij = P)
    }
    
    ## fi
  } else if( type == 2 ){
    
    param = matrix(0, ncol = 2, nrow = S)
    ## Attention unif not finite
    for (i in 1:S){
      
      ### ATTENTION UNIF not fait
      if ( distr[i] == "unif" ){
        ## censoring at the end
        if( cens.end == 1 && cens.beg == 0){
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          # vect.Nik = as.vector(t(Nik))
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*dunif(1:(Kmax), min = 0, max = par, log = TRUE))
               # + sum(Nei[i,]*log(1-sum(dpois(x = 0:(Kmax-1), lambda = par))))
               + sum(Nei[i,]*log(punif(q = 1:Kmax -1, min = 0, max = par, lower.tail = FALSE)))
               ## to faire pour les autres aussi !! 
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = Kmax, logvrais, method = "Brent", lower = 0, upper = 1e8)
          param[i,1] = CO2$par
          
          
          ## censoring at the beginning 
        }else if( cens.end == 0 && cens.beg == 1 ){   
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          vect.Nik = as.vector(t(Nik))
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*dunif(1:(Kmax), min = 0, max = par, log = TRUE))
               + sum(Nbi[i,]*log(punif(q = 1:Kmax -1, min = 0, max = par, lower.tail = FALSE)))
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = Kmax, logvrais, method = "Brent", lower = 0, upper = 1e8)
          param[i,1] = CO2$par
          
          ## censoring at the beginning and at the end
        }else if( cens.beg == 1 && cens.end == 1 ){
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*dunif(1:(Kmax), min = 0, max = par, log = TRUE))
               + sum(Nbi[i,]*log(punif(q = 1:Kmax -1, min = 0, max = par, lower.tail = FALSE)))
               + sum(Nei[i,]*log(punif(q = 1:Kmax -1, min = 0, max = par, lower.tail = FALSE)))
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = Kmax, logvrais, method = "Brent", lower = 0, upper = 1e8)
          param[i,1] = CO2$par
          
          ## No censoring  
        }else{        
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*dunif(1:(Kmax), min = 0, max = par, log = TRUE))
            )
          }
          
          ### Log-likelihood optimisation
          # CO2 = optim(par = 2, logvrais, method = "BFGS")
          # => done directely
          x = max.col(Nik, ties.method = "last")
          param[i,1] = x[i]
          
        }
        
      }else if ( distr[i] == "pois" ){
        ## censoring at the end
        if( cens.end == 1 && cens.beg == 0){
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          # vect.Nik = as.vector(t(Nik))
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*dpois(0:(Kmax-1), lambda = par, log = TRUE))
               # + sum(Nei[i,]*log(1-sum(dpois(x = 0:(Kmax-1), lambda = par))))
               + sum(Nei[i,]*log(ppois(q = 1:Kmax -1, par, lower.tail = FALSE)))
               ## to faire pour les autres aussi !! 
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = 10, logvrais, method = "Brent", lower = 0, upper = 1e8)
          param[i,1] = CO2$par
          
          
          ## censoring at the beginning 
        }else if( cens.end == 0 && cens.beg == 1 ){   
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          vect.Nik = as.vector(t(Nik))
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*dpois(0:(Kmax-1), lambda = par, log = TRUE))
               + sum(Nbi[i,]*log(ppois(q = 1:Kmax -1, par, lower.tail = FALSE)))
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = 3, logvrais, method = "Brent", lower = 0, upper = 1e8)
          param[i,1] = CO2$par
          
          ## censoring at the beginning and to the end
        }else if( cens.beg == 1 && cens.end == 1 ){
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*dpois(0:(Kmax-1), lambda = par, log = TRUE))
               + sum(Nbi[i,]*log(ppois(q = 1:Kmax -1, par, lower.tail = FALSE)))
               + sum(Nei[i,]*log(ppois(q = 1:Kmax -1, par, lower.tail = FALSE)))
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = 3, logvrais, method = "Brent", lower = 0, upper = 1e8)
          param[i,1] = CO2$par
          
          ## No censoring  
        }else{        
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*dpois(0:(Kmax-1), lambda = par, log = TRUE))
            )
          }
          
          ### Log-likelihood optimisation
          # CO2 = optim(par = 2, logvrais, method = "BFGS")
          # => faire directement
          x = rep(0:(Kmax-1), Nik[i,])
          param[i,1] = (sum(x)/sum(Nik[i,]))
          
        }
        
        
      }else if (distr[i] == "geom"){
        ## censoring at the end
        if( cens.end == 1 && cens.beg == 0){
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          vect.Nik = as.vector(t(Nik))
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*dgeom(0:(Kmax-1), prob = par, log = TRUE))
               + sum(Nei[i,]*log(pgeom(q = 1:Kmax -1, par, lower.tail = FALSE)))
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = 1/2, logvrais, method = "Brent", lower = 0, upper = 1)
          
          ### Get the parameters
          theta = CO2$par
          param[i,1] = theta
          
          ## censoring at the beginning 
        }else if( cens.end == 0 && cens.beg == 1 ){   
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*dgeom(0:(Kmax-1), prob = par, log = TRUE))
               + sum(Nbi[i,]*log(pgeom(q = 1:Kmax -1, par, lower.tail = FALSE)))
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = 1/2, logvrais, method = "Brent", lower = 0, upper = 1)
          
          ### Get the parameters
          theta = CO2$par
          param[i,1] = theta
          
          ## censoring at the beginning and to the end
        }else if( cens.beg == 1 && cens.end == 1 ){
          
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          vect.Nik = as.vector(t(Nik))
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*dgeom(0:(Kmax-1), prob = par, log = TRUE))
               + sum(Nbi[i,]*log(pgeom(q = 1:Kmax -1, par, lower.tail = FALSE)))
               + sum(Nei[i,]*log(pgeom(q = 1:Kmax -1, par, lower.tail = FALSE)))
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = 1/2, logvrais, method = "Brent", lower = 0, upper = 1)
          
          ### Get the parameters
          theta = CO2$par
          param[i,1] = theta
          
          ## No censoring  
        }else{        
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*dgeom(0:(Kmax-1), prob = par, log = TRUE))
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = 1/2, logvrais, method = "Brent", lower = 0, upper = 1)
          
          ### Get the parameters
          theta = CO2$par
          param[i,1] = theta
          
        }
        
      }else if ( distr[i] == "nbinom"){
        
        
        ## censoring at the end
        if( cens.end == 1 && cens.beg == 0){
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*dnbinom(x = 0:(Kmax-1), size = par[1], mu = par[2], log = TRUE))
               + sum(Nei[i,]*log(pnbinom(q = 1:Kmax -1, size = par[1], mu = par[2], lower.tail = FALSE)))
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = c(2, 2), logvrais, method = "Nelder-Mead")
          
          ### Get the parameters
          param[i,] = CO2$par
          
          
          ## censoring at the beginning 
        }else if( cens.end == 0 && cens.beg == 1 ){   
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*dnbinom(x = 0:(Kmax-1), size = par[1], mu = par[2], log = TRUE))
               + sum(Nbi[i,]*log(pnbinom(q = 1:Kmax -1, size = par[1], mu = par[2], lower.tail = FALSE)))
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = c(2,2), logvrais, method = "Nelder-Mead")
          
          ### Get the parameters
          param[i,] = CO2$par
          
          ## censoring at the beginning and at the end
        }else if( cens.beg == 1 && cens.end == 1 ){
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*dnbinom(x = 0:(Kmax-1), size = par[1], mu = par[2], log = TRUE))
               + sum(Nbi[i,]*log(pnbinom(q = 1:Kmax -1, size = par[1], mu = par[2], lower.tail = FALSE)))
               + sum(Nei[i,]*log(pnbinom(q = 1:Kmax -1, size = par[1], mu = par[2], lower.tail = FALSE)))
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = c(2,2), logvrais, method = "Nelder-Mead")
          
          ### Get the parameters
          param[i,] = CO2$par
          
          ## No censoring  
        }else{        
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          vect.Nik = as.vector(t(Nik))
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*dnbinom(x = 0:(Kmax-1), size = par[1], mu = par[2], log = TRUE))
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = c(2,2), logvrais, method = "Nelder-Mead")
          
          ### Get the parameters
          param[i,] = CO2$par
          
        }
        
      }else if ( distr[i] == "dweibull" ){
        
        ## censoring at the end
        if( cens.end == 1 && cens.beg == 0){
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*log(ddweibull(x = 1:Kmax, q = par[1], beta = par[2])))
               + sum(Nei[i,]*log(1-pdweibull(x = 1:Kmax, q = par[1], beta = par[2])))
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = c(0.5,1), logvrais, method = "Nelder-Mead")
          
          ### Get the parameters
          param[i,] = CO2$par
          
          ## censoring at the beginning 
        }else if( cens.end == 0 && cens.beg == 1 ){   
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*log(ddweibull(x = 1:Kmax, q = par[1], beta = par[2])))
               + sum(Nbi[i,]*log(1-pdweibull(x = 1:Kmax, q = par[1], beta = par[2])))
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = c(0.5,1), logvrais, method = "Nelder-Mead")
          
          ### Get the parameters
          param[i,] = CO2$par
          
          ## censoring at the beginning and to the end
        }else if( cens.beg == 1 && cens.end == 1 ){
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          vect.Nik = as.vector(t(Nik))
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*log(ddweibull(x = 1:Kmax, q = par[1], beta = par[2])))
               + sum(Nbi[i,]*log(1-sum(ddweibull(x = 1:Kmax, q = par[1], beta = par[2]))))
               + sum(Nei[i,]*log(1-pdweibull(x = 1:Kmax, q = par[1], beta = par[2])))
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = c(0.6,1), logvrais, method = "Nelder-Mead")
          
          ### Get the parameters
          param[i,] = CO2$par
          
          ## No censoring  
        }else{        
          
          
          #############################
          ## Estimation of the transition matrix
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Nik[i,]*log(ddweibull(x = 1:Kmax, q = par[1], beta = par[2])))
            )
          }
          
          ### Log-likelihood optimisation
          CO2 = optim(par = c(0.6, 1), logvrais, method = "Nelder-Mead")
          
          ### Get the parameters
          param[i,] = CO2$par 
        }
      }else{
        stop("The name of the distribution is not valid")
      }
      
    }
    
    
    ## Semi-Markov kernel
    kernel = .kernel_param_fi(distr = distr, param = param, Kmax = Kmax, pij = P, S = S)
    q = kernel$q
    ## Inital law
    if(nbSeq > 1){
      init = Nstart/sum(Nstart)
    }else{
      ## Computation of the limit distribution
      init = .limit.law(q = q, pij = P)
    }
    
    ## fj  
  }else if ( type == 3 ){
    
    param = matrix(0, nrow = S, ncol = 2)
    ## Attention unif pas fini
    for (i in 1:S){
      ### ATTENTION UNIF pas fait
      
      if ( distr[i] == "unif" ){
        ## censoring at the end
        if( cens.end == 1 && cens.beg == 0){
          if (S == 2){
            #############################
            ## Estimation of parameters 
            #############################
            ## vect.Nijk is the vector of the array Nijk
            vect.Nijk = as.vector(Nijk)
            ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
            ## estimation function of  p_{jt+1v} and theta_{jt+1v}
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            # Nij0 = Nij0[c(etat.cens),]
            
            
            ### Separation of parameters for likelihood estimation
            
            Nij1 = Nij.vect[1]
            Nij2 = Nij.vect[2]
            
            puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
            
           ## Functions to estimate
            logvrais = function(par){
              
              puv.vect0 = as.vector(t(puv))
              puv.vect = puv.vect0[-El.diag]
              # Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                 + sum(Njkv *dunif(1:(Kmax), min = 0, max = par[1], log = TRUE))
                 # + sum( as.vector(t(Nei))* log(1-(puv*sum(dpois(x = 0:(Kmax-1), lambda = par[(S*(S-2)+1):(S*(S-1))])))))
                 + sum( as.vector(t(Nei))* log(puv.vect*(punif(q = rep(1:Kmax -1,S), min = 0, max = par[1], lower.tail = FALSE))))
              )
            }
            
            ### Log-likelihood optimisation
            CO2<-optim(par=Kmax,logvrais, method ="Brent", lower = 0, upper = Kmax + 1)
            ## WARNINGS => production of NaN
            ### Get the parameters
            parC = CO2$par
            theta = parC[1]
            
            puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
            P = puv
            
            param[i,1] = theta
            
          }else{
            #############################
            ## Estimation of parameters 
            #############################
            ## Transformation of Nij pour manipuler les puv
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            ### Separation of parameters for likelihood estimation
            Nij1 = c()
            Nij2 = c()
            a = 1
            while (a < length(Nij.vect)){
              Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
              Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
              a = a + (S-1)
            }
            
            
           ## Functions to estimate
            logvrais = function(par){
              m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
              p = cbind(m, 1-rowSums(m))
              puv = rep(rowSums(p), each=Kmax)
              # Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                 + sum(Njkv *dunif(1:(Kmax), min = 0, max = par[(S*(S-2)+1)], log = TRUE))
                 # + sum( as.vector(t(Nei))* log(1-(puv*sum(dpois(x = 0:(Kmax-1), lambda = par[(S*(S-2)+1):(S*(S-1))])))))
                 + sum( as.vector(t(Nei))* log(puv*(punif(q = rep(1:Kmax -1,S), min = 0, max = par[(S*(S-2)+1)], lower.tail = FALSE))))
              )
            }
            
            ### CONSTRAINTS
            # puv > 0
            u0 = matrix(0, nrow = (S*(S-2)+1), ncol = (S*(S-2)+1))
            diag(u0) <- 1
            c0 = rep(0,(S*(S-2)+1))
            
            # puv < 1
            u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2)+1))
            diag(u1) <- -1
            c1 = c(rep(-1,S*(S-2)))
            
            u2 = rbind(u0,u1)
            c2 = c(c0,c1)
            
            ## verification of constraints
            # k2 <- length(c2)
            # n <- dim(u2)[2]
            # for(i in seq_len(k2)) {
            #   j <- which( u2[i,] != 0 )
            #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
            #   cat(" >= " )
            #   cat( c2[i], "\n" )
            # }
            
            ### Log-likelihood optimisation
            CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), Kmax),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
            ## WARNINGS => production of NaN
            ### Get the parameters
            parC = CO2$par
            parL = parC[1:(S*(S-2))]
            theta = parC[(S*(S-2)+1)]
            
            ### Transformation of parameters
            ## Add of parameters of type: pat = 1 - pac - pag
            ## in the vector par in optim
            b = 1
            parP = c()
            while( b <= length(parL) ){
              parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
              b = b +S-2
              
            }
            
            ## Add of zeros in diagonale for puv : 
            c = 1 
            parP0 = c()
            while( c < length(parP) ){
              parP0 = c(parP0, 0,parP[c:(c+(S-1))])
              c = c + S
            }
            parP0 = c(parP0, 0)
            P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
            P[which(is.na(P))] = 0
            param[i,1] = theta
            
          }
          
          ## censoring at the beginning 
        }else if( cens.end == 0 && cens.beg == 1 ){   
          
          #############################
          ## Estimation of puv
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          
          Nbjv = Nbj[i,]
          Njkv = Njk[i,]
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
         ## Functions to estimate
          logvrais = function(par){
            
            -( sum(Njkv *dunif(1:(Kmax), min = 0, max = par, log = TRUE))
               + sum( Nbjv * log(punif(q = 1:Kmax-1, min = 0, max = par, lower.tail = FALSE)))
            )
          }
          
          
          ### Log-likelihood optimisation
          CO2<-optim(par=Kmax,logvrais, method = "Brent", lower = 0, upper = 1e8)
          #### Warnings: NaN
          ### Get the parameters
          theta = CO2$par
          param[i,1] = theta
          
          
          ## censoring at the end and at the beginning
        }else if( cens.end == 1 && cens.beg == 1 ){
          
          if( S == 2 ){
            #############################
            ## Estimation of parameters 
            #############################
            ## vect.Nijk is the vector of the array Nijk
            vect.Nijk = as.vector(Nijk)
            ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
            ## estimation function of  p_{jt+1v} and theta_{jt+1v}
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            # Nij0 = Nij0[c(etat.cens),]
            
            
            ### Separation of parameters for likelihood estimation
            
            Nij1 = Nij.vect[1]
            Nij2 = Nij.vect[2]
            
            
           ## Functions to estimate
            logvrais = function(par){
              
              puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
              puv.vect0 = as.vector(t(puv))
              puv.vect = puv.vect0[-El.diag]
              # Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                 + sum(Njkv *dunif(1:(Kmax), min = 0, max = par[1], log = TRUE))
                 # + sum( as.vector(t(Nei))* log(1-(puv*sum(dpois(x = 0:(Kmax-1), lambda = par[(S*(S-2)+1):(S*(S-1))])))))
                 + sum( as.vector(t(Nei))* log(puv.vect*(punif(q = rep(1:Kmax -1,S), min = 0, max = par[1], lower.tail = FALSE))))
                 + sum( Nbjv * log(punif(q = 1:Kmax - 1, min = 0, max = par[1], lower.tail = FALSE)))
                 
              )
            }
            
            ### Log-likelihood optimisation
            CO2<-optim(par=Kmax,logvrais, method = "brent", lower = 0, upper = 1e9)
            ## WARNINGS => production of NaN
            ### Get the parameters
            parC = CO2$par
            theta = parC[1]
            param[i,1] = theta
            
            puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
            P = puv
            
          }else{
            #############################
            ## Estimation of parameters 
            #############################
            ## Transformation of Nij pour manipuler les puv
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            ### Separation of parameters for likelihood estimation
            ### l'estimation of puv dans la logvrais
            Nij1 = c()
            Nij2 = c()
            a = 1
            while (a < length(Nij.vect)){
              Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
              Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
              a = a + (S-1)
            }
            
            Nbjv = Nbj[i,]
            Njkv = Njk[i,] 
           ## Functions to estimate
            logvrais = function(par){
              m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
              p = cbind(m, 1-rowSums(m))
              puv = rep(rowSums(m), each = Kmax)
              Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                 + sum(Njkv *dunif(1:(Kmax), min = 0, max = par[(S*(S-2)+1)], log = TRUE))
                 + sum( as.vector(t(Nei))* log(puv*(punif(q = rep(1:Kmax -1,S), min = 0, max = par[(S*(S-2)+1)], lower.tail = FALSE))))
                 + sum( Nbjv * log(punif(q = 1:Kmax - 1, min = 0, max = par[(S*(S-2)+1)], lower.tail = FALSE)))
              )
            }
            
            ### CONSTRAINTS
            # puv > 0
            u0 = matrix(0, nrow = (S*(S-2)+1), ncol = (S*(S-2)+1))
            diag(u0) <- 1
            c0 = rep(0,(S*(S-2)+1))
            
            # puv < 1
            u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2)+1))
            diag(u1) <- -1
            c1 = c(rep(-1,S*(S-2)))
            
            u2 = rbind(u0,u1)
            c2 = c(c0,c1)
            
            ## verification of constraints
            # k2 <- length(c2)
            # n <- dim(u2)[2]
            # for(i in seq_len(k2)) {
            #   j <- which( u2[i,] != 0 )
            #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
            #   cat(" >= " )
            #   cat( c2[i], "\n" )
            # }
            
            ### Log-likelihood optimisation
            CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), Kmax),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
            #### Warnings : NaN
            ### Get the parameters
            parC = CO2$par
            parL = parC[1:(S*(S-2))]
            theta = parC[(S*(S-2)+1)]
            
            ### Transformation of parameters
            ##  Add of parameters of type: pat = 1 - pac - pag
            ## in the vector par in optim
            b = 1
            parP = c()
            while( b <= length(parL) ){
              parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
              b = b +S-2
              
            }
            
            ## Add of zeros in diagonale for puv : 
            c = 1 
            parP0 = c()
            while( c < length(parP) ){
              parP0 = c(parP0, 0,parP[c:(c+(S-1))])
              c = c + S
            }
            parP0 = c(parP0, 0)
            P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
            P[which(is.na(P))] = 0
            param[i,1] = theta
            
          }
          
          ## No censoring  
        }else{        
          
          #############################
          ## Estimation of puv
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          
          Nbjv = Nbj[i,]
          Njkv = Njk[i,]
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
         ## Functions to estimate
          logvrais = function(par){
            
            Nbjv = Nbj[i,]
            Njkv = Njk[i,]
            
            -( sum(Njkv *dunif(1:(Kmax), min = 0, max = par, log = TRUE))
            )
          }
          
          
          ### Log-likelihood optimisation
          # CO2<-optim(par=Kmax,logvrais, method = "Brent", lower = 0, upper = 1e8)
          #### Warnings : NaN
          ### Get the parameters
          # theta = CO2$par
          x = max.col(Njk, ties.method = "last")
          param[i,1] = x[i]
        }
      }else if ( distr[i] == "pois" ){
        
        ## censoring at the end
        if( cens.end == 1 && cens.beg == 0){
          if (S == 2){
            #############################
            ## Estimation of parameters 
            #############################
            ## vect.Nijk is the vector of the array Nijk
            vect.Nijk = as.vector(Nijk)
            ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
            ## estimation function of  p_{jt+1v} and theta_{jt+1v}
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            # Nij0 = Nij0[c(etat.cens),]
            
            
            ### Separation of parameters in order to compute log-likelihood
            Nij1 = Nij.vect[1]
            Nij2 = Nij.vect[2]
            
            puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
            
           ## Functions to estimate
            logvrais = function(par){
              
              puv.vect0 = as.vector(t(puv))
              puv.vect = puv.vect0[-El.diag]
              # Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                 + sum(Njkv *dpois(0:(Kmax-1), lambda = par[1], log = TRUE))
                 # + sum( as.vector(t(Nei))* log(1-(puv*sum(dpois(x = 0:(Kmax-1), lambda = par[(S*(S-2)+1):(S*(S-1))])))))
                 + sum( as.vector(t(Nei))* log(puv.vect*(ppois(q = rep(1:Kmax -1,S), par[1], lower.tail = FALSE))))
              )
            }
            
            ### Log-likelihood optimisation
            CO2<-optim(par=c(4),logvrais, method ="Brent", lower = 0, upper = Kmax)
            ## WARNINGS => production of NaN
            ### Get the parameters
            parC = CO2$par
            theta = parC[1]
            
            puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
            P = puv
            
            param[i,1] = theta
            
          }else{
            #############################
            ## Estimation of parameters 
            #############################
            ## Transformation of Nij pour manipuler les puv
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            ### Separation of parameters for likelihood estimation
            ### l'estimation of puv dans la logvrais
            Nij1 = c()
            Nij2 = c()
            a = 1
            while (a < length(Nij.vect)){
              Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
              Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
              a = a + (S-1)
            }
            
            
           ## Functions to estimate
            logvrais = function(par){
              m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
              p = cbind(m, 1-rowSums(m))
              puv = rep(rowSums(p), each=Kmax)
              # Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                 + sum(Njkv *dpois(0:(Kmax-1), lambda = par[(S*(S-2)+1)], log = TRUE))
                 # + sum( as.vector(t(Nei))* log(1-(puv*sum(dpois(x = 0:(Kmax-1), lambda = par[(S*(S-2)+1):(S*(S-1))])))))
                 + sum( as.vector(t(Nei))* log(puv*(ppois(q = rep(1:Kmax -1,S), par[(S*(S-2)+1)], lower.tail = FALSE))))
              )
            }
            
            ### CONSTRAINTS
            # puv > 0
            u0 = matrix(0, nrow = (S*(S-2)+1), ncol = (S*(S-2)+1))
            diag(u0) <- 1
            c0 = rep(0,(S*(S-2)+1))
            
            # puv < 1
            u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2)+1))
            diag(u1) <- -1
            c1 = c(rep(-1,S*(S-2)))
            
            u2 = rbind(u0,u1)
            c2 = c(c0,c1)
            
            ## verification of constraints
            # k2 <- length(c2)
            # n <- dim(u2)[2]
            # for(i in seq_len(k2)) {
            #   j <- which( u2[i,] != 0 )
            #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
            #   cat(" >= " )
            #   cat( c2[i], "\n" )
            # }
            
            ### Log-likelihood optimisation
            CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), 4),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
            ## WARNINGS => production of NaN
            ### Get the parameters
            parC = CO2$par
            parL = parC[1:(S*(S-2))]
            theta = parC[(S*(S-2)+1)]
            
            ### Transformation of parameters
            ## Add parameters of type  : pat = 1 - pac - pag
            ## in par vector of optim
            b = 1
            parP = c()
            while( b <= length(parL) ){
              parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
              b = b +S-2
              
            }
            
            ## Add zeros on diagonal for puv : 
            c = 1 
            parP0 = c()
            while( c < length(parP) ){
              parP0 = c(parP0, 0,parP[c:(c+(S-1))])
              c = c + S
            }
            parP0 = c(parP0, 0)
            P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
            P[which(is.na(P))] = 0
            param[i,1] = theta
            
          }
          
          ## censoring at the beginning 
        }else if( cens.end == 0 && cens.beg == 1 ){   
          
          #############################
          ## Estimation of puv
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          
          Nbjv = Nbj[i,]
          Njkv = Njk[i,]
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          ## Functions to estimate
          logvrais = function(par){

            -( sum(Njkv *dpois(0:(Kmax-1), lambda = par, log = TRUE))
               + sum( Nbjv * log(ppois(q = 1:Kmax-1, lambda = par, lower.tail = FALSE)))
            )
          }
          
          
          ### Log-likelihood optimisation
          CO2<-optim(par=5,logvrais, method = "Brent", lower = 0, upper = 1e8)
          #### Warnings : NaN
          ### Get the parameters
          theta = CO2$par
          param[i,1] = theta
          
          
          ## censoring at the end and at the beginning
        }else if( cens.end == 1 && cens.beg == 1 ){
          
          if( S == 2 ){
            #############################
            ## Estimation of parameters 
            #############################
            ## vect.Nijk is the vector of the array Nijk
            vect.Nijk = as.vector(Nijk)
            ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
            ## estimation function of  p_{jt+1v} and theta_{jt+1v}
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            # Nij0 = Nij0[c(etat.cens),]
            
            
            ### Separation of parameters for likelihood estimation
            
            Nij1 = Nij.vect[1]
            Nij2 = Nij.vect[2]
            
            
           ## Functions to estimate
            logvrais = function(par){
              
              puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
              puv.vect0 = as.vector(t(puv))
              puv.vect = puv.vect0[-El.diag]
              # Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                 + sum(Njkv *dpois(0:(Kmax-1), lambda = par[1], log = TRUE))
                 # + sum( as.vector(t(Nei))* log(1-(puv*sum(dpois(x = 0:(Kmax-1), lambda = par[(S*(S-2)+1):(S*(S-1))])))))
                 + sum( as.vector(t(Nei))* log(puv.vect*(ppois(q = rep(1:Kmax -1,S), par[1], lower.tail = FALSE))))
                 + sum( Nbjv * log(ppois(q = 1:Kmax - 1, lambda = par[1], lower.tail = FALSE)))
                 
              )
            }
            
            ### Log-likelihood optimisation
            CO2<-optim(par=c(4),logvrais, method = "brent", lower = 0, upper = 1e9)
            ## WARNINGS => production of NaN
            ### Get the parameters
            parC = CO2$par
            theta = parC[1]
            param[i,1] = theta
            
            puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
            P = puv
            
          }else{
            #############################
            ## Estimation of parameters 
            #############################
            ## Transformation of Nij for modifying puv
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            ### Separation of parameters in order to compute log-likelihood
            Nij1 = c()
            Nij2 = c()
            a = 1
            while (a < length(Nij.vect)){
              Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
              Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
              a = a + (S-1)
            }
            
            Nbjv = Nbj[i,]
            Njkv = Njk[i,] 
           ## Functions to estimate
            logvrais = function(par){
              m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
              p = cbind(m, 1-rowSums(m))
              puv = rep(rowSums(m), each = Kmax)
              Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                 + sum(Njkv *dpois(0:(Kmax-1), lambda = par[(S*(S-2)+1)], log = TRUE))
                 + sum( as.vector(t(Nei))* log(puv*(ppois(q = rep(1:Kmax -1,S), par[(S*(S-2)+1)], lower.tail = FALSE))))
                 + sum( Nbjv * log(ppois(q = 1:Kmax - 1, lambda = par[(S*(S-2)+1)], lower.tail = FALSE)))
              )
            }
            
            ### CONSTRAINTS
            # puv > 0
            u0 = matrix(0, nrow = (S*(S-2)+1), ncol = (S*(S-2)+1))
            diag(u0) <- 1
            c0 = rep(0,(S*(S-2)+1))
            
            # puv < 1
            u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2)+1))
            diag(u1) <- -1
            c1 = c(rep(-1,S*(S-2)))
            
            u2 = rbind(u0,u1)
            c2 = c(c0,c1)
            
            ## verification of constraints
            # k2 <- length(c2)
            # n <- dim(u2)[2]
            # for(i in seq_len(k2)) {
            #   j <- which( u2[i,] != 0 )
            #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
            #   cat(" >= " )
            #   cat( c2[i], "\n" )
            # }
            
            ### Log-likelihood optimisation
            CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), 4),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
            #### Warnings : NaN
            ### Get the parameters
            parC = CO2$par
            parL = parC[1:(S*(S-2))]
            theta = parC[(S*(S-2)+1)]
            
            ### Transformation of parameters
            ## Add parameters of type  : pat = 1 - pac - pag
            ## in par vector of optim
            b = 1
            parP = c()
            while( b <= length(parL) ){
              parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
              b = b +S-2
              
            }
            
            ## Add zeors on diagonal for puv : 
            c = 1 
            parP0 = c()
            while( c < length(parP) ){
              parP0 = c(parP0, 0,parP[c:(c+(S-1))])
              c = c + S
            }
            parP0 = c(parP0, 0)
            P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
            P[which(is.na(P))] = 0
            param[i,1] = theta
            
          }
          
          ## No censoring  
        }else{        
          
          #############################
          ## Estimation of puv
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          
          Nbjv = Nbj[i,]
          Njkv = Njk[i,]
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
         ## Functions to estimate
          logvrais = function(par){
            
            Nbjv = Nbj[i,]
            Njkv = Njk[i,]
            
            -( sum(Njkv *dpois(0:(Kmax-1), lambda = par, log = TRUE))
            )
          }
          
          
          ### Log-likelihood optimisation
          CO2<-optim(par=4,logvrais, method = "Brent", lower = 0, upper = 1e8)
          #### Warnings : NaN
          ### Get the parameters
          theta = CO2$par
          param[i,1] = theta
        }
        
      }else if ( distr[i] == "geom" ){
        
        ## censoring at the end
        if( cens.end == 1 && cens.beg == 0){
          if ( S == 2 ){
            #############################
            ## Estimation of parameters 
            #############################
            ## vect.Nijk is the vector of the array Nijk
            vect.Nijk = as.vector(Nijk)
            ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
            ## estimation function of  p_{jt+1v} and theta_{jt+1v}
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            # Nij0 = Nij0[c(etat.cens),]
            
            
            ### Separation of parameters for likelihood estimation
            
            Nij1 = Nij.vect[1]
            Nij2 = Nij.vect[2]
            
            
           ## Functions to estimate
            logvrais = function(par){
              
              puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
              puv.vect0 = as.vector(t(puv))
              puv.vect = puv.vect0[-El.diag]
              # Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                 + sum(Njkv *dgeom(0:(Kmax-1), prob = par[1], log = TRUE))
                 # + sum( as.vector(t(Nei))* log(1-(puv*sum(dpois(x = 0:(Kmax-1), lambda = par[(S*(S-2)+1):(S*(S-1))])))))
                 + sum( as.vector(t(Nei))* log(puv.vect*(pgeom(q = rep(1:Kmax -1,S), par[1], lower.tail = FALSE))))
              )
            }
            
            ### Log-likelihood optimisation
            CO2<-optim(par=c(0.4),logvrais, method ="Brent", lower = 0, upper = 1)
            ## WARNINGS => production of NaN
            ### Get the parameters
            parC = CO2$par
            theta = parC[1]
            
            puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
            P = puv
            
            param[i,1] = theta
          }else{
            #############################
            ## Estimation of parameters 
            #############################
            ## Transformation of Nij pour manipuler les puv
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            ### Separation of parameters in order to compute log-likelihood
            Nij1 = c()
            Nij2 = c()
            a = 1
            while (a < length(Nij.vect)){
              Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
              Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
              a = a + (S-1)
            }
            
            
           ## Functions to estimate
            logvrais = function(par){
              m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
              p = cbind(m, 1-rowSums(m))
              puv = rep(rowSums(m), each = Kmax)
              # Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                 + sum(Njkv *dgeom(0:(Kmax-1), prob = par[(S*(S-2)+1)], log = TRUE))
                 # + sum( as.vector(t(Nei))* log(1-(puv*sum(dgeom(x = 0:(Kmax-1), prob = par[(S*(S-2)+1):(S*(S-1))])))))
                 + sum( as.vector(t(Nei))* log(puv*(pgeom(q = rep(1:Kmax -1,S), par[(S*(S-2)+1)], lower.tail = FALSE))))
                 
              )
            }
            
            ### CONSTRAINTS
            # puv > 0
            u0 = matrix(0, nrow = (S*(S-2)+1), ncol = (S*(S-2)+1))
            diag(u0) <- 1
            c0 = rep(0,(S*(S-2)+1))
            
            # puv < 1
            u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2)+1))
            diag(u1) <- -1
            c1 = c(rep(-1,S*(S-2)))
            
            u2 = rbind(u0,u1)
            c2 = c(c0,c1)
            
            
            ## verification of constraints
            # k2 <- length(c2)
            # n <- dim(u2)[2]
            # for(i in seq_len(k2)) {
            #   j <- which( u2[i,] != 0 )
            #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
            #   cat(" >= " )
            #   cat( c2[i], "\n" )
            # }
            ### Log-likelihood optimisation
            CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), 1/4),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
            ## WARNINGS => production of NaN
            ### Get the parameters
            parC = CO2$par
            parL = parC[1:(S*(S-2))]
            theta = parC[(S*(S-2)+1)]
            
            ### Transformation of parameters
            ##  Add of parameters of type  : pat = 1 - pac - pag
            ## in the vector par in optim
            b = 1
            parP = c()
            while( b <= length(parL) ){
              parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
              b = b +S-2
              
            }
            
            ## Add of zeros in diagonale for puv : 
            c = 1 
            parP0 = c()
            while( c < length(parP) ){
              parP0 = c(parP0, 0,parP[c:(c+(S-1))])
              c = c + S
            }
            parP0 = c(parP0, 0)
            P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
            P[which(is.na(P))] = 0
            param[i,1] = theta
            
          }
          
          ## censoring at the beginning 
        }else if( cens.end == 0 && cens.beg == 1 ){   
          
          #############################
          ## Estimation of puv
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          Nbjv = Nbj[i,]
          Njkv = Njk[i,]
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
         ## Functions to estimate
          logvrais = function(par){
            
            Nbjv = Nbj[i,]
            Njkv = Njk[i,]
            
            -( sum(Njkv *dgeom(0:(Kmax-1), prob = par, log = TRUE))
               + sum( Nbjv * log(pgeom(q = 1:Kmax-1, prob = par, lower.tail = FALSE)))
            )
          }
          
          
          ### Log-likelihood optimisation
          CO2<-optim(par=1/4,logvrais, method = "Brent", lower = 0, upper = 1)
          #### Warnings : NaN
          ### Get the parameters
          theta = CO2$par
          param[i,1] = theta
          
          
          ## censoring at the end and at the beginning
        }else if( cens.end == 1 && cens.beg == 1 ){
          
          if ( S == 2 ){
            #############################
            ## Estimation of parameters 
            #############################
            ## vect.Nijk is the vector of the array Nijk
            vect.Nijk = as.vector(Nijk)
            ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
            ## estimation function of  p_{jt+1v} and theta_{jt+1v}
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            # Nij0 = Nij0[c(etat.cens),]
            
            
            ### Separation of parameters in order to compute log-likelihood 
            Nij1 = Nij.vect[1]
            Nij2 = Nij.vect[2]
            
            
           ## Functions to estimate
            logvrais = function(par){
            
              puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)  
              puv.vect0 = as.vector(t(puv))
              puv.vect = puv.vect0[-El.diag]
              # Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                 + sum(Njkv *dgeom(0:(Kmax-1), prob = par[1], log = TRUE))
                 # + sum( as.vector(t(Nei))* log(1-(puv*sum(dpois(x = 0:(Kmax-1), lambda = par[(S*(S-2)+1):(S*(S-1))])))))
                 + sum( as.vector(t(Nei))* log(puv.vect*(pgeom(q = rep(1:Kmax -1,S), par[1], lower.tail = FALSE))))
                 + sum( Nbjv * log(pgeom(q = 1:Kmax - 1, par[1], lower.tail = FALSE)))
                 
              )
            }
            
            ### Log-likelihood optimisation => change "brent"
            CO2<-optim(par=c(0.4),logvrais, method ="Brent", lower = 0, upper = 1)
            ## WARNINGS => production of NaN
            ### Get the parameters
            parC = CO2$par
            theta = parC[1]
            
            puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
            P = puv
            
            param[i,1] = theta
            
          }else{
            #############################
            ## Estimation of parameters 
            #############################
            ## Transformation of Nij pour manipuler les puv
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            ### Separation of parameters in order to compute log-likelihood
            Nij1 = c()
            Nij2 = c()
            a = 1
            while (a < length(Nij.vect)){
              Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
              Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
              a = a + (S-1)
            }
            
            Nbjv = Nbj[i,]
            Njkv = Njk[i,]
            
            
           ## Functions to estimate
            logvrais = function(par){
              m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
              p = cbind(m, 1-rowSums(m))
              puv = rep(rowSums(m), each = Kmax)
              Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                 + sum(Njkv *dgeom(0:(Kmax-1), prob = par[(S*(S-2)+1)], log = TRUE))
                 + sum( as.vector(t(Nei))* log(puv*(pgeom(q = rep(1:Kmax -1,S), par[(S*(S-2)+1)], lower.tail = FALSE))))
                 + sum( Nbjv * log(pgeom(q = 1:Kmax-1, prob = par[(S*(S-2)+1)], lower.tail = FALSE)))
              )
            }
            
            ### CONSTRAINTS
            # puv > 0
            u0 = matrix(0, nrow = (S*(S-2)+1), ncol = (S*(S-2)+1))
            diag(u0) <- 1
            c0 = rep(0,(S*(S-2)+1))
            
            # puv < 1
            u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2)+1))
            diag(u1) <- -1
            c1 = c(rep(-1,S*(S-2)))
            
            u2 = rbind(u0,u1)
            c2 = c(c0,c1)
            
            ## verification of constraints
            # k2 <- length(c2)
            # n <- dim(u2)[2]
            # for(i in seq_len(k2)) {
            #   j <- which( u2[i,] != 0 )
            #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
            #   cat(" >= " )
            #   cat( c2[i], "\n" )
            # }
            
            ### Log-likelihood optimisation
            CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), 1/3),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
            #### Warnings : NaN
            ### Get the parameters
            parC = CO2$par
            parL = parC[1:(S*(S-2))]
            theta = parC[(S*(S-2)+1)]
            
            ### Transformation of parameters
            ##  Add of parameters of type  : pat = 1 - pac - pag
            ## in the vector par in optim
            b = 1
            parP = c()
            while( b <= length(parL) ){
              parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
              b = b +S-2
              
            }
            
            ## Add zeros in diagonal for puv : 
            c = 1 
            parP0 = c()
            while( c < length(parP) ){
              parP0 = c(parP0, 0,parP[c:(c+(S-1))])
              c = c + S
            }
            parP0 = c(parP0, 0)
            P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
            P[which(is.na(P))] = 0
            param[i,1] = theta
            
          }
          
          ## No censoring  
        }else{        
          
          #############################
          ## Estimation of puv
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          Nbjv = Nbj[i,]
          Njkv = Njk[i,]
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
         ## Functions to estimate
          logvrais = function(par){
            
            Nbjv = Nbj[i,]
            Njkv = Njk[i,]
            
            -( sum(Njkv *dgeom(0:(Kmax-1), prob = par, log = TRUE))
            )
          }
          
          ### Log-likelihood optimisation
          CO2<-optim(par=1/4,logvrais, method = "Brent", lower = 0, upper = 1)
          #### Warnings : NaN
          ### Get the parameters
          theta = CO2$par
          param[i,1] = theta
        }
        
      }else if ( distr[i] == "nbinom"){
        
        ## censoring at the end
        if( cens.end == 1 && cens.beg == 0){
          
          if ( S == 2 ){
            #############################
            ## Estimation of parameters 
            #############################
            ## vect.Nijk is the vector of the array Nijk
            vect.Nijk = as.vector(Nijk)
            ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
            ## estimation function of  p_{jt+1v} and theta_{jt+1v}
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            # Nij0 = Nij0[c(etat.cens),]
            
            
            ### Separation of parameters for likelihood estimation
            
            Nij1 = Nij.vect[1]
            Nij2 = Nij.vect[2]
            
            
           ## Functions to estimate
            logvrais = function(par){
            
              puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
              puv.vect0 = as.vector(t(puv))
              puv.vect = puv.vect0[-El.diag]
              # Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                 + sum(Njkv *dnbinom(0:(Kmax-1), size = par[1], mu = par[2], log = TRUE))
                 # + sum( as.vector(t(Nei))* log(1-(puv*sum(dpois(x = 0:(Kmax-1), lambda = par[(S*(S-2)+1):(S*(S-1))])))))
                 + sum( as.vector(t(Nei))* log(puv.vect*(pnbinom(q = rep(1:Kmax -1,S), size = par[1], mu = par[2], lower.tail = FALSE))))
                 # + sum( Nbjv * log(pnbinom(q = 1:Kmax - 1, size = par[1], mu = par[2], lower.tail = FALSE)))
                 
              )
            }
            
            ### Log-likelihood optimisation => changer "brent"
            CO2<-optim(par=c(2,4),logvrais, method ="Nelder-Mead")
            ## WARNINGS => production of NaN
            ### Get the parameters
            parC = CO2$par
            alpha = parC[1]
            theta = parC[2]
            
            puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
            P = puv
            
            param[i,1] = alpha
            param[i,2] = theta
            
          }else{
            #############################
            ## Estimation of parameters 
            #############################
            ## Transformation of Nij pour manipuler les puv
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            ### Separation of parameters in order to compute log-likelihood
            Nij1 = c()
            Nij2 = c()
            a = 1
            while (a < length(Nij.vect)){
              Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
              Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
              a = a + (S-1)
            }
            
            
           ## Functions to estimate
            logvrais = function(par){
              m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
              p = cbind(m, 1-rowSums(m))
              puv = rep(rowSums(m), each = Kmax)
              # Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                 + sum(Njkv *dnbinom(0:(Kmax-1), size = par[(S*(S-2)+1)], mu = par[(S*(S-2)+2)], log = TRUE))
                 # + sum( as.vector(t(Nei))* log(1-(puv*sum(dnbinom(x = 0:(Kmax-1), size = par[(S*(S-2)+1):(S*(S-1)+1)], mu = par[(S*(S-1)+1):((S-1)*(S+1))])))))
                 + sum( as.vector(t(Nei))* log(puv*(pnbinom(q = rep(1:Kmax -1,S), size = par[(S*(S-2)+1)], mu = par[(S*(S-2)+2)], lower.tail = FALSE))))
              )
            }
            
            ### CONSTRAINTS
            # puv > 0
            u0 = matrix(0, nrow = (S*(S-2)+2), ncol = (S*(S-2)+2))
            diag(u0) <- 1
            c0 = rep(0,(S*(S-2)+2))
            
            # puv < 1
            u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2)+2))
            diag(u1) <- -1
            c1 = c(rep(-1,S*(S-2)))
            
            u2 = rbind(u0,u1)
            c2 = c(c0,c1)
            
            ## verification of constraints
            # k2 <- length(c2)
            # n <- dim(u2)[2]
            # for(i in seq_len(k2)) {
            #   j <- which( u2[i,] != 0 )
            #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
            #   cat(" >= " )
            #   cat( c2[i], "\n" )
            # }
            
            ### Log-likelihood optimisation
            CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), 2,4),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
            ## WARNINGS => production of NaN
            ### Get the parameters
            parC = CO2$par
            parL = parC[1:(S*(S-2))]
            alpha = parC[(S*(S-2)+1)]
            theta = parC[(S*(S-2)+2)]
            
            ### Transformation of parameters
            ##  Add of parameters of type  : pat = 1 - pac - pag
            ## in the vector par in optim
            b = 1
            parP = c()
            while( b <= length(parL) ){
              parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
              b = b +S-2
              
            }
            
            ## Add of zeros in diagonale for puv : 
            c = 1 
            parP0 = c()
            while( c < length(parP) ){
              parP0 = c(parP0, 0,parP[c:(c+(S-1))])
              c = c + S
            }
            parP0 = c(parP0, 0)
            P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
            P[which(is.na(P))] = 0
            param[i,1] = alpha
            param[i,2] = theta  
          }
          
          
          ## censoring at the beginning 
        }else if( cens.end == 0 && cens.beg == 1 ){   
          
          #############################
          ## Estimation of puv
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          Nbjv = Nbj[i,]
          Njkv = Njk[i,]
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
         ## Functions to estimate
          logvrais = function(par){
            
            Nbjv = Nbj[i,]
            Njkv = Njk[i,]
            
            -( sum(Njkv *dnbinom(0:(Kmax-1), size = par[1], mu = par[2], log = TRUE))
               + sum( Nbjv * log(pnbinom(q = 1:Kmax-1, size = par[1], mu = par[2], lower.tail = FALSE)))
            )
          }
          
          
          ### Log-likelihood optimisation
          CO2<-optim(par=c(2,4),logvrais)
          #### Warnings : NaN
          ### Get the parameters
          param[i,] = CO2$par
          
          ## censoring at the end and at the beginning
        }else if( cens.end == 1 && cens.beg == 1 ){
          
          if ( S == 2 ){
            #############################
            ## Estimation of parameters 
            #############################
            ## vect.Nijk is the vector of the array Nijk
            vect.Nijk = as.vector(Nijk)
            ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
            ## estimation function of  p_{jt+1v} and theta_{jt+1v}
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            # Nij0 = Nij0[c(etat.cens),]
            
            
            ### Separation of parameters for likelihood estimation
            
            Nij1 = Nij.vect[1]
            Nij2 = Nij.vect[2]
            
            
           ## Functions to estimate
            logvrais = function(par){
              
              puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
              puv.vect0 = as.vector(t(puv))
              puv.vect = puv.vect0[-El.diag]
              # Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                 + sum(Njkv *dnbinom(0:(Kmax-1), size = par[1], mu = par[2], log = TRUE))
                 + sum( as.vector(t(Nei))* log(puv*(pnbinom(q = rep(1:Kmax -1,S), size = par[1], mu = par[2], lower.tail = FALSE))))
                 # + sum( as.vector(t(Nei))* log(1-(puv*sum(dnbinom(x = 0:(Kmax-1), size = par[(S*(S-2)+1):(S*(S-1)+1)], mu = par[(S*(S-1)+1):((S-1)*(S+1))])))))
                 + sum( Nbjv * log(pnbinom(q = 1:Kmax-1, size = par[1], mu = par[2], lower.tail = FALSE)))
              )
              
            }
            
            ### Log-likelihood optimisation => change "brent"
            CO2<-optim(par=c(2,4),logvrais, method ="Nelder-Mead")
            ## WARNINGS => production of NaN
            ### Get the parameters
            parC = CO2$par
            alpha = parC[1]
            theta = parC[2]
            
            puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
            P = puv
            
            param[i,1] = alpha
            param[i,2] = theta
            
          }else{
            #############################
            ## Estimation of parameters 
            #############################
            ## Transformation of Nij pour manipuler les puv
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            ### Separation of parameters in order to compute log-likelihood
            Nij1 = c()
            Nij2 = c()
            a = 1
            while (a < length(Nij.vect)){
              Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
              Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
              a = a + (S-1)
            }
            
            Nbjv = Nbj[i,]
            Njkv = Njk[i,]
            
            
           ## Functions to estimate
            logvrais = function(par){
              m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
              p = cbind(m, 1-rowSums(m))
              puv = rep(rowSums(p), each = Kmax)
              Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                 + sum(Njkv *dnbinom(0:(Kmax-1), size = par[(S*(S-2)+1)], mu = par[(S*(S-2)+2)], log = TRUE))
                 + sum( as.vector(t(Nei))* log(puv*(pnbinom(q = rep(1:Kmax -1,S), size = par[(S*(S-2)+1)], mu = par[(S*(S-2)+2)], lower.tail = FALSE))))
                 # + sum( as.vector(t(Nei))* log(1-(puv*sum(dnbinom(x = 0:(Kmax-1), size = par[(S*(S-2)+1):(S*(S-1)+1)], mu = par[(S*(S-1)+1):((S-1)*(S+1))])))))
                 + sum( Nbjv * log(pnbinom(q = 1:Kmax-1, size = par[(S*(S-2)+1)], mu = par[(S*(S-2)+2)], lower.tail = FALSE)))
              )
              
            }
            
            ### CONSTRAINTS
            # puv > 0
            u0 = matrix(0, nrow = (S*(S-2)+2), ncol = (S*(S-2)+2))
            diag(u0) <- 1
            c0 = rep(0,(S*(S-2)+2))
            
            # puv < 1
            u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2)+2))
            diag(u1) <- -1
            c1 = c(rep(-1,S*(S-2)))
            
            u2 = rbind(u0,u1)
            c2 = c(c0,c1)
            
            ## verification of constraints
            # k2 <- length(c2)
            # n <- dim(u2)[2]
            # for(i in seq_len(k2)) {
            #   j <- which( u2[i,] != 0 )
            #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
            #   cat(" >= " )
            #   cat( c2[i], "\n" )
            # }
            
            ### Log-likelihood optimisation
            CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), 10,4),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
            ## WARNINGS => production of NaN
            ### Get the parameters
            parC = CO2$par
            parL = parC[1:(S*(S-2))]
            alpha = parC[(S*(S-2)+1)]
            theta = parC[(S*(S-2)+2)]
            
            ### Transformation of parameters
            ##  Add of parameters of type  : pat = 1 - pac - pag
            ## in the vector par in optim
            b = 1
            parP = c()
            while( b <= length(parL) ){
              parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
              b = b +S-2
              
            }
            
            ## Add of zeros in diagonale for puv : 
            c = 1 
            parP0 = c()
            while( c < length(parP) ){
              parP0 = c(parP0, 0,parP[c:(c+(S-1))])
              c = c + S
            }
            parP0 = c(parP0, 0)
            P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
            P[which(is.na(P))] = 0
            param[i,1] = alpha
            param[i,2] = theta
            
          }
          
          ## No censoring  
        }else{        
          
          #############################
          ## Estimation of puv
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          Nbjv = Nbj[i,]
          Njkv = Njk[i,]
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          ## Functions to estimate
          logvrais = function(par){
            
            Nbjv = Nbj[i,]
            Njkv = Njk[i,]
            
            -( sum(Njkv *dnbinom(0:(Kmax-1), size = par[1], mu = par[2], log = TRUE))
            )
          }
          
          
          ### Log-likelihood optimisation
          CO2<-optim(par=c(2,4),logvrais)
          #### Warnings : NaN
          ### Get the parameters
          param[i,] = CO2$par
        }
        
      }else if ( distr[i] == "dweibull" ){
        
        ## censoring at the end
        if( cens.end == 1 && cens.beg == 0){
          if ( S == 2 ){
            
            #############################
            ## Estimation of parameters 
            #############################
            ## vect.Nijk is the vector of the array Nijk
            vect.Nijk = as.vector(Nijk)
            ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
            ## estimation function of  p_{jt+1v} and theta_{jt+1v}
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            # Nij0 = Nij0[c(etat.cens),]
            
            
            ### Separation of parameters for likelihood estimation
            
            Nij1 = Nij.vect[1]
            Nij2 = Nij.vect[2]
            
            puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
            
           ## Functions to estimate
            logvrais = function(par){
              
              puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
              puv.vect0 = as.vector(t(puv))
              puv.vect = puv.vect0[-El.diag]
              # Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                 + sum(Njkv *log(ddweibull(1:Kmax, q = par[1], beta = par[2])))
                 # + sum( as.vector(t(Nei))* log(1-(puv*sum(dpois(x = 0:(Kmax-1), lambda = par[(S*(S-2)+1):(S*(S-1))])))))
                 + sum( as.vector(t(Nei))* log(puv*(1-pdweibull(x = rep(1:Kmax,S), q = par[(S*(S-2)+1)], beta = par[(S*(S-2)+2)]))))
                 # + sum( Nbjv * log(pnbinom(q = 1:Kmax - 1, size = par[1], mu = par[2], lower.tail = FALSE)))
                 
              )
            }
            
            ### Log-likelihood optimisation => change "brent"
            CO2<-optim(par=c(0.6,0.8),logvrais, method ="Nelder-Mead")
            ## WARNINGS => production of NaN
            ### Get the parameters
            parC = CO2$par
            q = parC[1]
            beta = parC[2]
            
            puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
            P = puv
            
            param[i,1] = q
            param[i,2] = beta
            
            
          }else{
            #############################
            ## Estimation of parameters 
            #############################
            ## Transformation of Nij pour manipuler les puv
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            ### Separation of parameters in order to compute log-likelihood
            Nij1 = c()
            Nij2 = c()
            a = 1
            while (a < length(Nij.vect)){
              Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
              Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
              a = a + (S-1)
            }
            
            
            ## Functions to estimate
            logvrais = function(par){
              m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
              p = cbind(m, 1-rowSums(m))
              puv = rep(rowSums(p), each=Kmax)
              # Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                 + sum(Njkv*log(ddweibull(1:(Kmax), q = par[(S*(S-2)+1)], beta = par[(S*(S-2)+2)])))
                 + sum( as.vector(t(Nei))* log(puv*(1-pdweibull(x = rep(1:Kmax,S), q = par[(S*(S-2)+1)], beta = par[(S*(S-2)+2)]))))
                 # + sum( as.vector(t(Nei))* log(1-(puv*sum(ddweibull(x = 1:(Kmax), q = par[(S*(S-2)+1):(S*(S-1)+1)], beta = par[(S*(S-1)+1):((S-1)*(S+1))])))))
              )
            }
            
            
            ### CONSTRAINTS
            # puv > 0
            u0 = matrix(0, nrow = (S*(S-2)+2), ncol = (S*(S-2)+2))
            diag(u0) <- 1
            c0 = rep(0,(S*(S-2)+2))
            
            # puv < 1
            u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2)+2))
            diag(u1) <- -1
            c1 = c(rep(-1,S*(S-2)))
            
            u2 = rbind(u0,u1)
            c2 = c(c0,c1)
            
            ## verification of constraints
            # k2 <- length(c2)
            # n <- dim(u2)[2]
            # for(i in seq_len(k2)) {
            #   j <- which( u2[i,] != 0 )
            #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
            #   cat(" >= " )
            #   cat( c2[i], "\n" )
            # }
            
            ### Log-likelihood optimisation
            CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), 0.6,0.8),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
            ## WARNINGS => production of NaN
            ### Get the parameters
            parC = CO2$par
            parL = parC[1:(S*(S-2))]
            q = parC[(S*(S-2)+1)]
            beta = parC[(S*(S-2)+2)]
            
            ### Transformation of parameters
            ##  Add of parameters of type  : pat = 1 - pac - pag
            ## in the vector par in optim
            b = 1
            parP = c()
            while( b <= length(parL) ){
              parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
              b = b +S-2
              
            }
            
            ## Add of zeros in diagonale for puv : 
            c = 1 
            parP0 = c()
            while( c < length(parP) ){
              parP0 = c(parP0, 0,parP[c:(c+(S-1))])
              c = c + S
            }
            parP0 = c(parP0, 0)
            P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
            P[which(is.na(P))] = 0
            param[i,1] = q
            param[i,2] = beta
            
          }
          
          ## censoring at the beginning 
        }else if( cens.end == 0 && cens.beg == 1 ){   
          
          #############################
          ## Estimation of puv
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          Nbjv = Nbj[i,]
          Njkv = Njk[i,]
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
         ## Functions to estimate
          logvrais = function(par){
            
            Nbjv = Nbj[i,]
            Njkv = Njk[i,]
            
            -( sum(Njkv *log(ddweibull(1:(Kmax), q = par[1], beta = par[2])))
               + sum( Nbjv * log(1-pdweibull(x = 1:Kmax, q = par[1], beta = par[2])))
            )
          }
          
          
          ### Log-likelihood optimisation
          CO2<-optim(par=c(0.6,0.8),logvrais)
          #### Warnings : NaN
          ### Get the parameters
          param[i,] = CO2$par
          
          ## censoring at the end and at the beginning
        }else if( cens.end == 1 && cens.beg == 1 ){
          
          if ( S == 2 ){
            #############################
            ## Estimation of parameters 
            #############################
            ## vect.Nijk is the vector of the array Nijk
            vect.Nijk = as.vector(Nijk)
            ##  Nijk[i,j,] correspond for Nijk vector to vect.Nijk[i+S*(j-1)+plus]
            ## estimation function of  p_{jt+1v} and theta_{jt+1v}
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            # Nij0 = Nij0[c(etat.cens),]
            
            
            ### Separation of parameters for likelihood estimation
            
            Nij1 = Nij.vect[1]
            Nij2 = Nij.vect[2]
            
            
           ## Functions to estimate
            logvrais = function(par){
              
              puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
              puv.vect0 = as.vector(t(puv))
              puv.vect = puv.vect0[-El.diag]
              # Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
     
              -( sum(Nij.vect*log(puv.vect)) #sum(par[1:(S-2)])))
                 + sum(Njkv *log(ddweibull(1:(Kmax), q = par[1], beta = par[2])))
                 # + sum( as.vector(t(Nei))* log(1-puv*sum(ddweibull(x = 1:(Kmax), q = par[(S*(S-2)+1):(S*(S-1)+1)], beta = par[(S*(S-1)+1):((S-1)*(S+1))]))))
                 + sum( as.vector(t(Nei))* log(puv*(1-pdweibull(x = rep(1:Kmax,S), q = par[1], beta = par[2]))))
                 + sum( Nbjv * log(1-pdweibull(x = 1:(Kmax), q = par[1], beta = par[2])))
              )
              
            }
            
            ### Log-likelihood optimisation => change "brent"
            CO2<-optim(par=c(0.6,0.8),logvrais, method ="Nelder-Mead")
            ## WARNINGS => production of NaN
            ### Get the parameters
            parC = CO2$par
            q = parC[1]
            beta = parC[2]
            
            puv = matrix(c(0,1,1,0), nrow = S, byrow = TRUE)
            P = puv
            
            param[i,1] = q
            param[i,2] = beta
            
          }else{
            #############################
            ## Estimation of parameters 
            #############################
            ## Transformation of Nij pour manipuler les puv
            Nij.vect0<-as.vector(t(Nij)) #[as.vector(t(Nij))!=0]
            El.diag = seq(from= 1, to = S*S, by = S+1)
            Nij.vect = Nij.vect0[-El.diag]
            ### Separation of parameters in order to compute log-likelihood
            Nij1 = c()
            Nij2 = c()
            a = 1
            while (a < length(Nij.vect)){
              Nij1 = c(Nij1,Nij.vect[a:(a+(S-3))]) ## correspond to the (S-2) first values of Nij on each line
              Nij2 = c(Nij2, Nij.vect[a+(S-2)]) ## correspond to the last value of Nij on each line
              a = a + (S-1)
            }
            
            Nbjv = Nbj[i,]
            Njkv = Njk[i,]
            
            
            ## Functions to estimate
            logvrais = function(par){
              m = matrix(par[1:(S*(S-2))], ncol = S-2, byrow = TRUE)
              puv = rep(rowSums(p), each=Kmax)
              Nbjv = Nbj[i,]
              Njkv = Njk[i,]
              
              -( sum(Nij1*log(par[1:(S*(S-2))])) + sum(Nij2*log(1-rowSums(m))) #sum(par[1:(S-2)])))
                 + sum(Njkv *log(ddweibull(1:(Kmax), q = par[(S*(S-2)+1)], beta = par[(S*(S-2)+2)])))
                 # + sum( as.vector(t(Nei))* log(1-puv*sum(ddweibull(x = 1:(Kmax), q = par[(S*(S-2)+1):(S*(S-1)+1)], beta = par[(S*(S-1)+1):((S-1)*(S+1))]))))
                 + sum( as.vector(t(Nei))* log(puv*(1-pdweibull(x = rep(1:Kmax,S), q = par[(S*(S-2)+1)], beta = par[(S*(S-2)+2)]))))
                 + sum( Nbjv * log(1-pdweibull(x = 1:(Kmax), q = par[(S*(S-2)+1)], beta = par[(S*(S-2)+2)])))
              )
              
            }
            
            
            ### CONSTRAINTS
            # puv > 0
            u0 = matrix(0, nrow = (S*(S-2)+2), ncol = (S*(S-2)+2))
            diag(u0) <- 1
            c0 = rep(0,(S*(S-2)+2))
            
            # puv < 1
            u1 = matrix(0, nrow = S*(S-2), ncol = (S*(S-2)+2))
            diag(u1) <- -1
            c1 = c(rep(-1,S*(S-2)))
            
            u2 = rbind(u0,u1)
            c2 = c(c0,c1)
            
            ## verification of constraints
            # k2 <- length(c2)
            # n <- dim(u2)[2]
            # for(i in seq_len(k2)) {
            #   j <- which( u2[i,] != 0 )
            #   cat(paste( u2[i,j], " * ", "x[", (1:n)[j], "]", sep="", collapse=" + " ))
            #   cat(" >= " )
            #   cat( c2[i], "\n" )
            # }
            
            ### Log-likelihood optimisation
            CO2<-constrOptim(theta=c(rep(1/(S-1),S*(S-2)), 0.6,0.8),logvrais, ui=u2,ci=c2, method ="Nelder-Mead")
            ## error : initial parameters !!!
            ### Get the parameters
            parC = CO2$par
            parL = parC[1:(S*(S-2))]
            q = parC[(S*(S-2)+1)]
            beta = parC[(S*(S-2)+2)]
            
            ### Transformation of parameters
            ##  Add of parameters of type  : pat = 1 - pac - pag
            ## in the vector par in optim
            b = 1
            parP = c()
            while( b <= length(parL) ){
              parP = c(parP, parL[b:(b+(S-3))], 1-sum(parL[b:(b+(S-3))]))
              b = b +S-2
              
            }
            
            ## Add of zeros on diagonal for puv : 
            c = 1 
            parP0 = c()
            while( c < length(parP) ){
              parP0 = c(parP0, 0,parP[c:(c+(S-1))])
              c = c + S
            }
            parP0 = c(parP0, 0)
            P = matrix(parP0, nrow = S, ncol = S, byrow = TRUE)
            P[which(is.na(P))] = 0
            param[i,1] = q
            param[i,2] = beta
            
            
          }
          
          ## No censoring  
        }else{        
          
          #############################
          ## Estimation of puv
          #############################
          P = Nij/Ni
          P[which(is.na(P))] = 0
          
          Nbjv = Nbj[i,]
          Njkv = Njk[i,]
          
          #############################
          ## Estimation of parameters  of the distribution 
          #############################
          ## Functions to estimate
          logvrais = function(par){
            
            Nbjv = Nbj[i,]
            Njkv = Njk[i,]
            
            -( sum(Njkv *log(ddweibull(1:(Kmax), q = par[1], beta = par[2])))
            )
          }
          
          
          ### Log-likelihood optimisation
          CO2<-optim(par=c(0.6,0.8),logvrais)
          #### Warnings : NaN
          ### Get the parameters
          param[i,] = CO2$par
        } 
        
      }else{
        stop("The name of the distribution is not valid")
      }
      
    }
    
    ## Semi-Markov kernel
    kernel = .kernel_param_fj(distr = distr, param = param, Kmax = Kmax, pij = P, S = S)
    q = kernel$q
    ## Inital law
    if(nbSeq > 1){
      init = Nstart/sum(Nstart)
    }else{
      ## calcul of the limit distribution
      init = .limit.law(q = q, pij = P)
    }
    
    ## f
  }else if( type == 4) {
    
    param = vector(mode = "numeric", length = 2)
    
    ### ATTENTION UNIF pas fait
    
    if ( distr == "unif" ){
      ## censoring at the end
      if( cens.end == 1 && cens.beg == 0){
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*dunif(1:(Kmax), min = 0, max = par, log = TRUE))
             + sum(Ne*log(punif(q = 1:Kmax-1, min = 0, max = par, lower.tail = FALSE)))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = Kmax, logvrais, method = "Brent", lower = 0, upper = Kmax + 1)
        
        
        ### Get the parameters
        theta = CO2$par
        param[1] = theta
        
        
        ## censoring at the beginning 
      }else if( cens.end == 0 && cens.beg == 1 ){   
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*dunif(1:(Kmax), min = 0, max = par, log = TRUE))
             + sum(Nb*log(punif(q = 1:Kmax-1, min = 0, max = par, lower.tail = FALSE)))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = Kmax, logvrais, method = "Brent", lower = 0, upper = 1e8)
        
        
        ### Get the parameters
        theta = CO2$par
        param[1] = theta
        
        ## censoring at the beginning and to the end
      }else if( cens.beg == 1 && cens.end == 1 ){
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*dunif(1:(Kmax), min = 0, max = par, log = TRUE))
             + sum(Neb*log(punif(q = 1:Kmax-1, min = 0, max = par, lower.tail = FALSE)))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = Kmax, logvrais, method = "Brent", lower = 0, upper = 1e8)
        
        ### Get the parameters
        theta = CO2$par
        param[1] = theta
        
        ## No censoring  
      }else{        
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*dunif(1:(Kmax), min = 0, max = par, log = TRUE))
          )
        }
        
        ### Log-likelihood optimisation
        # CO2 = optim(par = 2, logvrais, method = "Brent", lower = 0, upper = 10000)
        # => faire directement
        x = Nk[!is.null(Nk)]
        theta = length(x)
        
        ### Get the parameters
        # theta = CO2$par
        param[1] = theta 
      }
      
      
    }else if ( distr == "pois" ){
      ## censoring at the end
      if( cens.end == 1 && cens.beg == 0){
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*dpois(0:(Kmax-1), lambda = par, log = TRUE))
             + sum(Ne*log(ppois(q = 1:Kmax-1, lambda = par, lower.tail = FALSE)))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = 10, logvrais, method = "Brent", lower = 0, upper = 10000)
        
        
        ### Get the parameters
        theta = CO2$par
        param[1] = theta
        
        
        ## censoring at the beginning 
      }else if( cens.end == 0 && cens.beg == 1 ){   
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*dpois(0:(Kmax-1), lambda = par, log = TRUE))
             + sum(Nb*log(ppois(q = 1:Kmax-1, lambda = par, lower.tail = FALSE)))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = 3, logvrais, method = "Brent", lower = 0, upper = 1e8)
        
        
        ### Get the parameters
        theta = CO2$par
        param[1] = theta
        
        ## censoring at the beginning and to the end
      }else if( cens.beg == 1 && cens.end == 1 ){
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*dpois(0:(Kmax-1), lambda = par, log = TRUE))
             + sum(Neb*log(ppois(q = 1:Kmax-1, lambda = par, lower.tail = FALSE)))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = 3, logvrais, method = "Brent", lower = 0, upper = 1e8)
        
        ### Get the parameters
        theta = CO2$par
        param[1] = theta
        
        ## No censoring  
      }else{        
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*dpois(0:(Kmax-1), lambda = par, log = TRUE))
          )
        }
        
        ### Log-likelihood optimisation
        # CO2 = optim(par = 2, logvrais, method = "Brent", lower = 0, upper = 10000)
        # => faire directement
        x = rep(0:(Kmax-1), Nk)
        theta = (sum(x)/sum(Nk))
        
        ### Get the parameters
        # theta = CO2$par
        param[1] = theta 
      }
      
      
    }else if (distr == "geom"){
      ## censoring at the end
      if( cens.end == 1 && cens.beg == 0){
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
       ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*dgeom(0:(Kmax-1), prob = par, log = TRUE))
             + sum(Ne*log(pgeom(q = 1:Kmax-1, prob = par, lower.tail = FALSE)))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = 1/(S-1), logvrais, method = "Brent", lower = 0, upper = 1)
        
        ### Get the parameters
        theta = CO2$par
        param[1] = theta
        
        ## censoring at the beginning 
      }else if( cens.end == 0 && cens.beg == 1 ){   
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*dgeom(0:(Kmax-1), prob = par, log = TRUE))
             + sum(Nb*log(pgeom(q = 1:Kmax-1, prob = par, lower.tail = FALSE)))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = 1/2, logvrais, method = "Brent", lower = 0, upper = 1)
        
        ### Get the parameters
        theta = CO2$par
        param[1] = theta
        
        ## censoring at the beginning and to the end
      }else if( cens.beg == 1 && cens.end == 1 ){
        
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*dgeom(0:(Kmax-1), prob = par, log = TRUE))
             + sum(Neb*log(pgeom(q = 1:Kmax-1, prob = par, lower.tail = FALSE)))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = 1/2, logvrais, method = "Brent", lower = 0, upper = 1)
        
        ### Get the parameters
        theta = CO2$par
        param[1] = theta
        
        ## No censoring  
      }else{        
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*dgeom(0:(Kmax-1), prob = par, log = TRUE))
             # + sum(Neb*log(1-sum(dgeom(x = 0:(Kmax-1), prob = par))))
          )
        }
        
        ### Log-likelihood optimisation
        # CO2 = optim(par = 1/(S-1), logvrais, method = "Brent", lower = 0, upper = 1)
        x = rep(0:(Kmax-1), Nk)
        theta = (sum(x)/sum(Nk))
        
        ### Get the parameters
        # theta = CO2$par
        param[1] = 1/theta
        
      }
      
    }else if ( distr == "nbinom"){
      
      
      ## censoring at the end
      if( cens.end == 1 && cens.beg == 0){
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*dnbinom(x = 0:(Kmax-1), size = par[1], mu = par[2], log = TRUE))
             + sum(Ne*log(pnbinom(q = 1:Kmax-1, size = par[1], mu = par[2], lower.tail = FALSE)))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = c(2, 2), logvrais, method = "Nelder-Mead")
        
        ### Get the parameters
        param = CO2$par
        
        
        ## censoring at the beginning 
      }else if( cens.end == 0 && cens.beg == 1 ){   
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*dnbinom(x = 0:(Kmax-1), size = par[1], mu = par[2], log = TRUE))
             + sum(Nb*log(pnbinom(q = 1:Kmax-1, size = par[1], mu = par[2], lower.tail = FALSE)))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = c(2, 2), logvrais, method = "Nelder-Mead")
        
        ### Get the parameters
        param = CO2$par
        
        
        ## censoring at the beginning and at the end
      }else if( cens.beg == 1 && cens.end == 1 ){
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*dnbinom(x = 0:(Kmax-1), size = par[1], mu = par[2], log = TRUE))
             + sum(Neb*log(pnbinom(q = 1:Kmax-1, size = par[1], mu = par[2], lower.tail = FALSE)))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = c(2, 2), logvrais, method = "Nelder-Mead")
        
        ### Get the parameters
        param = CO2$par
        
        
        ## No censoring  
      }else{        
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*dnbinom(x = 0:(Kmax-1), size = par[1], mu = par[2], log = TRUE))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = c(2, 2), logvrais, method = "Nelder-Mead")
        
        ### Get the parameters
        param = CO2$par
        
      }
      
    }else if ( distr == "dweibull" ){
      
      ## censoring at the end
      if( cens.end == 1 && cens.beg == 0){
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*log(ddweibull(x = 1:Kmax, q = par[1], beta = par[2])))
             + sum(Ne*log(1-pdweibull(x = 1:Kmax, q = par[1], beta = par[2])))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = c(0.5,1), logvrais, method = "Nelder-Mead")
        
        ### Get the parameters
        param = CO2$par
        
        ## censoring at the beginning 
      }else if( cens.end == 0 && cens.beg == 1 ){   
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*log(ddweibull(x = 1:Kmax, q = par[1], beta = par[2])))
             + sum(Nb*log(1-pdweibull(x = 1:Kmax, q = par[1], beta = par[2])))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = c(0.5,1), logvrais, method = "Nelder-Mead")
        
        ### Get the parameters
        param = CO2$par
        
        ## censoring at the beginning and to the end
      }else if( cens.beg == 1 && cens.end == 1 ){
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*log(ddweibull(x = 1:Kmax, q = par[1], beta = par[2])))
             + sum(Neb*log(1-pdweibull(x = 1:Kmax, q = par[1], beta = par[2])))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = c(0.5,1), logvrais, method = "Nelder-Mead")
        
        ### Get the parameters
        param = CO2$par
        
        ## No censoring  
      }else{        
        
        
        #############################
        ## Estimation of the transition matrix
        #############################
        P = Nij/Ni
        P[which(is.na(P))] = 0
        
        #############################
        ## Estimation of parameters  of the distribution 
        #############################
        
        ## Functions to estimate
        logvrais = function(par){
          
          -( sum(Nk*log(ddweibull(x = 1:Kmax, q = par[1], beta = par[2])))
          )
        }
        
        ### Log-likelihood optimisation
        CO2 = optim(par = c(0.5,1), logvrais, method = "Nelder-Mead")
        
        ### Get the parameters
        param = CO2$par
      }
    }else{
      stop("The name of the distribution is not valid")
    }
    
    ## semi-markovian kernel
    kernel = .kernel_param_f(distr = distr, param = param, Kmax = Kmax, pij = P, S = S)
    q = kernel$q
    ## Inital law
    if(nbSeq > 1){
      init = Nstart/sum(Nstart)
    }else{
      ## calcul of the distribution  limite
      init = .limit.law(q = q, pij = P)
    }
    
  }else{
    stop("The parameter \"TypeSojournTime\" must be only \"fij\", \"fi\", \"fj\" or \"f\".")
  }
  
  ## Question : which parameters to return ? 
  .write.results(TypeSojournTime, cens.beg, cens.end, P, param)
  
  q[which(is.na(q))]<-0
  
  return(list(init = init, Ptrans = P, param = param, q = q))
  
}
