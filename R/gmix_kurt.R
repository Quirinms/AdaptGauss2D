gmix_kurt = function(Parm, Data, Threshold=1.0, Debug=0){
  # V = gmix_kurt(Parm,Data,Threshold,Debug)
  # V$Parm
  #
  # DESCRIPTION
  # Subroutine to update Gaussian mixture (5 Operations):
  # 4. Splitting modes (gmix_kurt.R)
  #
  # See the matlab documentation for more information.
  #
  # Kurtosis-based mode splitting based on
  # N. Vlassis A. Likas,  "The Kurtosis-EM algorithm for
  # Gaussian mixture modelling, to be published IEEE SMC
  # Transactions, in 1999
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
  # Data[1:n,1:d]    Numerical matrix with normalized data. N samples with DIM
  #                  feature dimensions.
  #
  # OPTIONAL
  # Threshold         Numerical value. The threshold to determine if a mode
  #                   should be split.
  # Debug             do some outputs. Default=0.
  #
  # OUTPUT
  # Parm              Updated Input parameters, see input description above.
  #
  #
  # in /dbt/EMforGauss/
  # Author QS 2021
  MAXMODE = 100
  N       = dim(Data)[1]
  DIM     = dim(Data)[2]
  wts     = Parm$Modes$Weight
  nmode   = length(wts)
  numColsCholesky = dim(Parm$Modes$Covariance)[2]
  # compute mode functions at each sample
  px = matrix(0, N, nmode)
  for(k in 1:nmode){
    idxFrom = numColsCholesky*(k-1)+1
    idxTo   = numColsCholesky*k
    tmpvar = Parm$Modes$Covariance[idxFrom:idxTo, ]
    V      = lqr_eval(Data, Parm$Modes$Mean[k,], tmpvar)
    px[,k] = V$LogPDF
  }
  # compute membership probabilities of modes for each sample
  if(nmode == 1){
    w = matrix(1, 1, N)
  }else{
    w = matrix(0, nmode, N)
    for(i in 1:nmode){
      w[i,] = (log(wts[i]) + px[,i])
    }
    mx = apply(w,2,max)
    for(i in 1:nmode){
      w[i,] = exp(w[i,] - mx)
    }
    nor = colSums(w)
    for(i in 1:nmode){
      w[i,] =  w[i,] / nor
    }
  }
  kall = c()
  jall = matrix(0, nmode, 1)
  mall = matrix(0, nmode, 1)
  for(k in 1:nmode){
    idxFrom = numColsCholesky*(k-1)+1
    idxTo   = numColsCholesky*k
    tmpvar  = Parm$Modes$Covariance[idxFrom:idxTo, ]
    # Get axes
    res_svd = svd(tmpvar)
    V       = res_svd$v
    S       = res_svd$d
    xt = Data - pracma::repmat(Parm$Modes$Mean[k,], N, 1)    # Remove the mean
    # for each dimension, project onto this vector
    mmax= -1
    for(j in 1:DIM){
      tmp = xt %*% V[,j]/S[j]
      mw  = max(w[k,])
      iw  = which(w[k,] > mw/2)
      sw = sum(w[k,])
      K1 = sum(w[k,] * tmp**1) / sw
      K2 = sum(w[k,] * tmp**2) / sw
      K3 = sum(w[k,] * tmp**3) / sw
      K4 = sum(w[k,] * tmp**4) / sw
      if(sw/30/DIM > 10){
        fac = 1
      }else if(sw/30/DIM > 2){
        fac = .75
      }else if(sw/30/DIM > 1){
        fac = .5
      }else{
        fac = 0
      }
      merit = 2*(abs(K4-3)  + abs(K3)) * fac * K2
      if(merit > Threshold){
        kall = c(kall, k)
      }
      if(Debug){
        print(paste(k,": ",j,", K2= ",K2,", K3= ",K3,", K4= ",K4,
                    ", sw= ",sw,", N = ",N,", merit = ",merit,
                    ", s = ",S[j],", fac = ",fac))
        #hist(tmp(iw)/S(j,j),64)
      }
      if(merit > mmax){
        mmax = merit
        jm = j
      }
    }
    mall[k] = mmax
    jall[k] = jm
  }
  ntot = nmode
  for(k in 1:nmode){
    if(mall[k] > Threshold){
      jm = jall[k]
      if(ntot >= MAXMODE){
        print("Kurt.m: Warning - MAXMODE exceeded")
        return()
      }
      if(Debug > 0){
        print(paste("KURT: Splitting a mode, dimension ", jm, " of ", DIM))
      }
      idxFrom = numColsCholesky*(k-1)+1
      idxTo   = numColsCholesky*k
      tmpvar = Parm$Modes$Covariance[idxFrom:idxTo, ]
      # Get Axes
      res_svd = svd(tmpvar)
      S       = res_svd$d
      V       = res_svd$v
      m1 = Parm$Modes$Mean[k,] + V[,jm] * S[jm]/2
      m2 = Parm$Modes$Mean[k,] - V[,jm] * S[jm]/2
      Parm$Modes$Weight[k] = wts[k]/2
      Parm$Modes$Weight    = c(Parm$Modes$Weight, wts[k]/2)
      Parm$Modes$Mean[k,]  = m1
      Parm$Modes$Mean      = rbind(Parm$Modes$Mean, m2)
      Parm$Modes$Covariance = rbind(Parm$Modes$Covariance, tmpvar)
      ntot = ntot + 1
      wts  = Parm$Modes$Weight
    }
  }
  rownames(Parm$Modes$Mean) = NULL
  return(list(Parm=Parm))
}
#
#
#
#
#
