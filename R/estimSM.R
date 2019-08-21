estimSM = function(file = NULL, seq = NULL, E, TypeSojournTime = "fij", distr = "NP", cens.end = 0, cens.beg = 0){
  #require(seqinr)
  #require(DiscreteWeibull)

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

  if (length(distr) == 1){
    if(distr == "NP"){
      if(cens.beg == 0 && cens.end == 0){
        .estimNP_aucuneCensure(file = file, seq, E, TypeSojournTime = TypeSojournTime, distr = "NP")
      }else if(cens.beg == 1 && cens.end == 0){
        .estimNP_CensureFin(file = file, seq, E, TypeSojournTime = TypeSojournTime)
      }else if(cens.beg == 0 && cens.end == 1){

      }else{

      }
    }else{
      .estim.plusTraj(seq = seq, E, TypeSojournTime = TypeSojournTime, distr = distr, cens.end = cens.end, cens.beg = cens.beg)
    }
  }else{ 
    .estim.plusTraj(seq = seq, E, TypeSojournTime = TypeSojournTime, distr = distr, cens.end = cens.end, cens.beg = cens.beg)
  }

}
