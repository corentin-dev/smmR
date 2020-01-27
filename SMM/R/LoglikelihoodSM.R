### INPUTS
## seq: list of sequences
## alphabet: vector of state space 
## mu: vector of initial distribution
## puv: matrix of trasition probabilities 
## fij: parameter of the distributions 
## distr: matrix of disributions
## TypeSojournTime: characer "fij", "fi", "fj" or "f"
LoglikelihoodSM = function(seq, E, mu, Ptrans, distr = "NP", param = NULL, laws = NULL, TypeSojournTime){

  lois = as.vector(distr)
  ## length of the state space
  S = length(E)
  if(dim(Ptrans)[1] != S || dim(Ptrans)[2] != S){
    stop("The size of the matrix Ptrans must be equal to SxS with S = length(E)")  
  }

  ## parametric case
  if(is.matrix(distr) || is.array(distr) || distr != "NP"){
    nonparametric=FALSE
    if(length(distr) == 1 && distr != "NP" && !is.null(laws)){
      stop("The parameter \"param\" must be used with the parameter \"distr\"")
    }
    if(is.null(param)){
    }
  }
  ## non parametric case
  else{
    nonparametric=TRUE
    if( nonparametric && is.null(laws) ){
      stop("The parameter \"param\" must be used with the parameter \"distr\"")
    }
  }

  if( is.null(param) && is.null(laws) ){
    stop("One of the two parameters \"param\" (with distr) or \"laws\" (without distr) should be given")
  }

  if(!is.list(seq)){
    stop("The parameter \"seq\" should be a list")
  }

  ## All sequences 
  nbSeq<-length(seq)
  Kmax = 0
  J<-list()
  T<-list()
  L<-list()
  M<-list()
  seqStart = vector(mode = "numeric", length = S)

  for (k in 1:nbSeq){
    ## Manipulations
    # seqStart[alphabet==seq[[k]][1]] = mu[alphabet==seq[[k]][1]] + 1 
    JT<-.donneJT(unlist(seq[[k]]))
    J[[k]]<-.code(JT[1,],E)
    T[[k]]<-as.numeric(JT[2,])
    L[[k]]<-T[[k]] - c(1,T[[k]][-length(T[[k]]-1)]) ## sojourn time
    Kmax <- max(Kmax, max(L[[k]])) 
  }

  ## J: list
  ## L: list
  ## S: alphabet size
  ## Kmax: maximal sojourn time
  ################
  ## Get all the counts
  res = .comptage(J = J, L = L, S = S, Kmax = Kmax)
  Nuv = res$LNij
  Nuvk = res$LNijk
  Nuk = res$LNik
  Nvk = res$LNjk
  Nk = res$LNk
  Nstarti = res$Nstarti
  Nstart = res$NStart

  LV = list()
  for (k in 1:nbSeq) {
    ## Comptuation of log-likelihood terms
    diago = seq(from = 1, to = S*S, by = S+1)
    Nuv.vect = Nuv[[k]][-diago]
    puv.vect = Ptrans[-diago]
    i = 1
    j = 0
    diag = rep(diago, Kmax)
    while (i < Kmax){
      diag[i:(i+(S-1))] = diag[i:(i+(S-1))] + j*S*S
      i = i + S
      j = j + 1
    }
    Nuvk.vect = Nuvk[[k]][-diag]

    if( "NP" %in% lois ){
      fij = laws
      if(TypeSojournTime == "fij"){
        fuv.vect = fij[-diag]
        Nuvk0 = Nuvk.vect[-which(is.infinite(fuv.vect))]
        fuv0 = fuv.vect[-which(is.infinite(Nuvk.vect))]
        LV[[k]] = sum(Nstart[[k]]*log(mu)) + sum(Nuv.vect*log(puv.vect)) + sum(Nuvk0*log(fuv0))

      }else if(TypeSojournTime == "fi"){
        ## Get only the first non-null element of each matrix
        FTEST = apply(X = fij, MARGIN = 2, FUN = unique)
        FT = apply(FTEST, 1, unique)
        fuv.vect = FT[-which(FT == 0)]      
        Nuvk0 = Nuvk.vect[-which(is.infinite(fuv.vect))]
        fuv0 = fuv.vect[-which(is.infinite(Nuvk.vect))]
        LV[[k]] = sum(Nstart[[k]]*log(mu)) + sum(Nuv.vect*log(puv.vect)) + sum(Nuvk0*log(fuv0))

      }else if(TypeSojournTime == "fj"){
        FTEST = apply(X = fij, MARGIN = 1, FUN = unique)
        FT = apply(FTEST, 1, unique)
        fuv.vect = FT[-which(FT == 0)]
        Nuvk0 = Nuvk.vect[-which(is.infinite(fuv.vect))]
        fuv0 = fuv.vect[-which(is.infinite(Nuvk.vect))]
        LV[[k]] = sum(Nstart[[k]]*log(mu)) + sum(Nuv.vect*log(puv.vect)) + sum(Nuvk0*log(fuv0))

      }else if(TypeSojournTime == "f"){
        FTEST = apply(X = fij, MARGIN = 3, FUN = unique)
        FT = apply(FTEST, 2, unique)
        fuv.vect = FT[-which(FT == 0)]
        Nuvk0 = Nuvk.vect[-which(is.infinite(fuv.vect))]
        fuv0 = fuv.vect[-which(is.infinite(Nuvk.vect))]
        LV[[k]] = sum(Nstart[[k]]*log(mu)) + sum(Nuv.vect*log(puv.vect)) + sum(Nuvk0*log(fuv0))

      }else{
        stop("TypeSojournTime must be equal to \"fij\", \"fi\", \"fj\" or \"f\" ")
      }

    }else{

      if(TypeSojournTime == "fij"){
        qfij = .kernel_param_fij(distr = distr, param = param, Kmax = Kmax, pij = Ptrans, S = S)
        fuv = qfij$f
        fuv.vect = fuv[-diag]
        Nuvk.vect = Nuvk.vect[-which(is.infinite(log(fuv.vect)))]
        fuv.vect = fuv.vect[-which(is.infinite(log(fuv.vect)))] 
        ## Log-likelihood
        LV[[k]] = sum(Nstart[[k]]*log(mu)) + sum(Nuv.vect*log(puv.vect)) + sum(Nuvk.vect*log(fuv.vect))

      }else if(TypeSojournTime == "fi"){
        qfi.j = .kernel_param_fi(distr = distr, param = param, Kmax = Kmax, pij = Ptrans, S = S)
        fu.v = qfij$f
        fuv.vect = fu.v[-diag]
        Nuvk.vect = Nuvk.vect[-which(is.infinite(log(fuv.vect)))]
        fuv.vect = fuv.vect[-which(is.infinite(log(fuv.vect)))] 
        ## Log-likelihood
        LV[[k]] = sum(Nstart[[k]]*log(mu)) + sum(Nuv.vect*log(puv.vect)) + sum(Nuvk.vect*log(fuv.vect))

      }else if( TypeSojournTime == "fj" ){
        qfi.j = .kernel_param_fj(distr = distr, param = param, Kmax = Kmax, pij = Ptrans, S = S)
        fu.v = qfij$f
        fuv.vect = fu.v[-diag]
        Nuvk.vect = Nuvk.vect[-which(is.infinite(log(fuv.vect)))]
        fuv.vect = fuv.vect[-which(is.infinite(log(fuv.vect)))] 
        ## Log-likelihood
        LV[[k]] = sum(Nstart[[k]]*log(mu)) + sum(Nuv.vect*log(puv.vect)) + sum(Nuvk.vect*log(fuv.vect))


      }else if( TypeSojournTime == "f" ){
        qf = .kernel_param_f(distr = distr, param = param, Kmax = Kmax, pij = Ptrans, S = S)
        fu.v = qf$f
        fuv.vect = fu.v[-diag]
        Nuvk.vect = Nuvk.vect[-which(is.infinite(log(fuv.vect)))]
        fuv.vect = fuv.vect[-which(is.infinite(log(fuv.vect)))] 
        ## Log-likelihood
        LV[[k]] = sum(Nstart[[k]]*log(mu)) + sum(Nuv.vect*log(puv.vect)) + sum(Nuvk.vect*log(fuv.vect))
      }else{
        stop("TypeSojournTime must be equal to \"fij\", \"fi\", \"fj\" or \"f\" ")
      }

    }  
  }


  return (list(L = LV, Kmax = Kmax))
}
