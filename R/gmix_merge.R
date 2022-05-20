gmix_merge = function(Parm, MaxCloseness, Verbose=0){
  # V = gmix_merge(Parm,MaxCloseness, Verbose=0)
  # V$parm
  # 
  # DESCRIPTION
  # Subroutine to update Gaussian mixture (5 Operations):
  # 3. Merging modes (gmix_merge.R)
  # This method creates a single mode from two nearly identical modes.
  # The closeness of two modes is determined by mode_dist.R.
  # See the matlab documentation for more information.
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
  # MaxCloseness     Numerical value: Maximum mode closeness. Use more negative
  #                  values to promote mode consolidation.  Use higher values
  #                  for larger dimension (suggest -0.5 times DIM).
  # 
  # OPTIONAL
  # Verbose          do some outputs , Default ==0;
  # 
  # OUTPUT
  # Parm             Updated Input parameters, see input description above.
  # 
  # 
  # in /dbt/EMforGauss/
  # Author QS 2021
  wts   = Parm$Modes$Weight
  nmode = length(wts)
  numColsCholesky = dim(Parm$Modes$Covariance)[2]
  if(nmode<2){
    return(list(Parm=Parm))
  }
  npurge = 0
  mx     = -1.0e100
  for(k1 in 1:nmode){
    if((k1+1) <= nmode){
      for(k2 in (k1+1):nmode){
        
        idxFrom1 = numColsCholesky*(k1-1)+1
        idxFrom2 = numColsCholesky*(k2-1)+1
        idxTo1   = numColsCholesky*k1
        idxTo2   = numColsCholesky*k2
        
        if((wts[k1]>0)&(wts[k2]>0)){
          V    = mode_dist(Parm$Modes$Mean[k1,],
                           Parm$Modes$Mean[k2,],
                           Parm$Modes$Covariance[idxFrom1:idxTo1,],
                           Parm$Modes$Covariance[idxFrom2:idxTo2,])
          dist = V$Dist
          if(dist > mx){
            kmax = c(k1, k2)
            mx   = dist
          }
          if((dist > MaxCloseness) & (npurge < ((nmode/2)+1))){
            npurge = npurge+1
            V = my_merge(wts[k1],
                        Parm$Modes$Mean[k1,],
                        Parm$Modes$Covariance[idxFrom1:idxTo1,],
                        wts[k2],
                        Parm$Modes$Mean[k2,],
                        Parm$Modes$Covariance[idxFrom2:idxTo2,])
            w1 = V$Weight
            m1 = V$Mean
            v1 = V$Covariance
            if(wts[k1] > wts[k2]){
              if(Verbose==1){
                print(paste("Merging ", k2, " into ", k1, ", d = ", dist, sep=""))
              }
              wts[k1] = w1
              Parm$Modes$Covariance[idxFrom1:idxTo1,] = v1
              Parm$Modes$Mean[k1,] = m1
              wts[k2] = 0
            }else{
              if(Verbose==1){
                print(paste("Merging ", k1, " into ", k2, ", d = ", dist, sep=""))
              }
              wts[k2] = w1
              Parm$Modes$Covariance[idxFrom2:idxTo2,] = v1
              Parm$Modes$Mean[k2,] = m1
              wts[k1] = 0
            }
          }
        }
        if(wts[k1] == 0){
          break
        }
      }
    }
  }
  for(i in 1:length(wts)){
    Parm$Modes$Weight[i] = wts[i]
  }
  if(Verbose == 1){
    print(paste("gmix_merge: CLOSEST are ", kmax[1], " and ", kmax[2], ", DIST = ", mx, sep=""))
  }
  V    = gmix_deflate(Parm, 1.0e-10, 1.0e-10, Verbose)
  Parm = V$Parm
  return(list(Parm=Parm))
}
#
#
#
#
#