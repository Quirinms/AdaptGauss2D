gmix_deflate = function(Parm, MinWeight1, MinWeightAll, Verbose=0){
  # Parm = gmix_deflate(Parm,MinWeight1,MinWeightAll)
  # 
  # DESCRIPTION
  # Subroutine to update Gaussian mixture (5 Operations):
  # 2. Pruning modes (gmix_deflate.R)
  # The method gmix_deflate.R is killing weak modes (a mode is another name for one
  # of the L mixture components). A weak mode is found by testing the correct mode
  # entry of weight vector wts to see if it falls below a threshold.
  # See the matlab documentation for more information.
  # Only one mode per run can be killed.
  # 
  # INPUT
  # Parm             Nested list with parameters for GMM.
  #                  Features carrying permanent values:
  #                  Features$Names [1:d]      String vector with feature names.
  #                  Features$MinStd [1:l]     Vector of covariance constraints.
  #                  Modes list of modifyable GMM parameters for all modes.
  #                  Modes$Covariance [1:d*l, 1:d]     Numerical matrix with l
  #                        many square matrices stacked vertically with the
  #                        covariance matrix.
  #                  Modes$Mean [1:l, 1:d]             Numerical matrix with l
  #                        row vectors containing the d dimensional means.
  #                  Modes$Weight [1:l]                Numerical vector Modes
  #                        weights for each mean.
  # MinWeight1       Numerical value: First mode below this threshold is killed.
  # MinWeightAll     Numerical value: All modes below this threshold are killed.
  # 
  # OPTIONAL
  # Verbose          do some outputs , Default=0;
  # 
  # OUTPUT
  # Parm             updated Input parameters, see above.
  # 
  # 
  # in /dbt/EMforGauss/
  # Author QS 2021
  wts   = Parm$Modes$Weight
  means = Parm$Modes$Mean
  nmode = length(wts)
  numColsCholesky = dim(Parm$Modes$Covariance)[2]
  if(nmode==1){
    return(list("Parm"=Parm))
  }
  ndefl=0
  # --- deflate weakest mode ---
  flag = 1
  while(flag == 1){
    ymn = min(wts)
    imn = which.min(wts)
    if((ymn < MinWeightAll) | ((ymn < MinWeight1) & (ndefl < 1))){
      if(ymn >= MinWeightAll){
        ndefl = ndefl+1
      }
      flag = 1
      if(Verbose == 1){
        print(paste("gmix_deflate: Deflating mode ", imn, ", weight = ", ymn, ",
                    nmode = ", nmode-1, sep=""))
      }
      ii = 1
      for(i in 1:nmode){
        if(i != imn){
          idxFromi  = numColsCholesky*(i-1)+1
          idxFromii = numColsCholesky*(ii-1)+1
          idxToi    = numColsCholesky*i
          idxToii   = numColsCholesky*ii
          Parm$Modes$Covariance[idxFromii:idxToii,] = Parm$Modes$Covariance[idxFromi:idxToi,]
          Parm$Modes$Mean[ii,]  = Parm$Modes$Mean[i,]
          Parm$Modes$Weight[ii] = Parm$Modes$Weight[i]
          ii = ii + 1
        }
      }
      nmode = nmode-1
      idxTo   = numColsCholesky*nmode
      Parm$Modes$Covariance = Parm$Modes$Covariance[1:idxTo,]
      Parm$Modes$Mean       = Parm$Modes$Mean[1:nmode,]
      Parm$Modes$Weight     = Parm$Modes$Weight[1:nmode]
      wts = Parm$Modes$Weight
      wts = wts/sum(wts)
    }else{
      flag = 0
    }
  }
  wts = Parm$Modes$Weight
  wts = wts/sum(wts)
  for(i in 1:length(wts)){
    Parm$Modes$Weight[i] = wts[i]
  }
  return(list(Parm=Parm))
}
#
#
#
#
#
