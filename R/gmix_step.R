gmix_step = function(Parm, Data, Bias=0, DataWTS=NULL){
  # V = gmix_step(Parm, Data, Bias=0, DataWTS=NULL)
  # V$Parm
  # V$Q
  #
  # DESCRIPTION
  # Subroutine to update Gaussian mixture (5 Operations):
  # 1. E-M algorithm (gmix_step.R)
  # E-M algorithm (expectation maximization algorithm)
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
  # Data[1:n,1:d]    Numerical matrix with normalized data. N samples with DIM
  #                  feature dimensions.
  #
  # OPTIONAL
  # Bias             Binary value: Covariance constraint method.
  #                  Choose: 1=BIAS, 0=CONSTRAINT
  # DataWTS[1:N]     Numerical Vector, which allows individually weighting input
  #                  data. Default: DataWTS = matrix(1,N,1)
  #
  # OUTPUT
  # Parm             Updated Input parameters, see input description above.
  # Q                Numerical value: Total log-likelihood output
  #                  (weighted by DataWTS).
  #
  #
  # Dr. Paul M. Baggenstoss
  # Naval Undersea Warfare Center
  # Newport, RI
  # p.m.baggenstoss@ieee.org
  #
  # in /dbt/EMforGauss/
  # Author QS 2021
  if(is.null(Parm)){
    warning("Input parameter Parm is not given. Returning.")
    return()
  }
  if(is.null(Data)){
    warning("Input parameter x is not given. Returning.")
    return()
  }
  N   = dim(Data)[1]
  DIM = dim(Data)[2]
  #nmode = length(Parm$modes)
  nmode = length(Parm$Modes$Weight)
  numColsCholesky = dim(Parm$Modes$Covariance)[2]
  if(is.null(DataWTS)){
    DataWTS = matrix(1, N, 1)
  }
  # compute mode functions at each sample
  px = matrix(0, N, nmode)
  for(k in 1:nmode){
    #mean = Parm[[paste("modes", k, sep="")]]$mean
    #cholesky_covar = Parm[[paste("modes", k, sep="")]]$cholesky_covar
    mean           = Parm$Modes$Mean[k,]
    idxFrom        = numColsCholesky*(k-1)+1
    idxTo          = numColsCholesky*k
    cholesky_covar = Parm$Modes$Covariance[idxFrom:idxTo,]
    V              = lqr_eval(Data, mean, cholesky_covar) # Verbraucht viel zeit!
    px[,k]         = V$LogPDF
  }
  if(nmode==1){
    Q = sum(px * DataWTS)
    w = DataWTS
    a = sum(w)
  }else{
    w = matrix(0, N, nmode)
    for(i in 1:nmode){
      weight = Parm$Modes$Weight[i]
      w[,i] = (log(weight) + px[,i])
    }
    Q = sum(log(rowSums(exp(w))) * DataWTS) # Verbraucht viel zeit
    mx = max(w)
    for(i in 1:nmode){
      w[,i] = exp( w[,i] - mx ) # Verbraucht viel zeit
    }
    nor = rowSums(w)
    eps = 10**(-5)
    nor[which(nor <= 0)] = eps # nor can have zero entries and thus must be corrected with eps
    for(i in 1:nmode){
      w[,i] = (DataWTS * w[,i]) / nor
    }
    a = colSums(w)
    a[which(a <= 0)] = eps 
  }
  # Jetzt gibts Q, w, a
  # Update means
  for(i in 1:nmode){
    Parm$Modes$Mean[i,] = (1/a[i]) * (t(Data) %*% (w[,i]))
  }
  # Update variance
  for(k in 1:nmode){
    uu = sqrt((w[,k]) / a[k])
    #mean = Parm[[paste("modes", k, sep="")]]$mean
    mean = Parm$Modes$Mean[k,]
    # Subtract mean and scale by weights
    tmpidx = (Data - pracma::repmat(mean, N, 1)) * pracma::repmat(as.matrix(uu), 1, DIM) # Verbraucht viel zeit
    #tmpidx[is.na(tmpidx)] = 0
    #print(tmpidx)
    #print(any(is.na(tmpidx)))
    res_qr = qr(tmpidx)
    q = qr.Q(res_qr)
    tmpvar = qr.R(res_qr)
    MinStd = Parm$Features$MinStd[1:2]
    if(Bias == 0){
      # constrains on the variances
      # first determine eigen decomp of C=R' * R = V * S^2 * V'
      res_svd = svd(tmpvar)
      S = res_svd$d
      U = res_svd$u
      V = res_svd$v
      S = pmax(S, sqrt(diag(t(V) %*% diag(MinStd^2) %*% V )))
      tmpvar = U %*% diag(S) %*% t(V)
      # here you must make it upper triangular
      res_qr = qr(tmpvar)
      q      = qr.Q(res_qr)
      tmpvar = qr.R(res_qr)
    }else{
      tmpidx = rbind(tmpvar, diag(MinStd))
      res_qr = qr(tmpidx)
      q      = qr.Q(res_qr)
      tmpvar = qr.R(res_qr)
    }
    #Parm[[paste("modes", k, sep="")]]$cholesky_covar = tmpvar
    #Parm[[paste("modes", k, sep="")]]$weight = a[k]
    idxFrom  = numColsCholesky*(k-1)+1
    idxTo    = numColsCholesky*k
    Parm$Modes$Covariance[idxFrom:idxTo,] = tmpvar
    Parm$Modes$Weight[k] = a[k]
  }
  wts = Parm$Modes$Weight
  wts = wts/sum(wts)
  wts[which(wts <= 0)] = 0
  for(i in 1:nmode){
    Parm$Modes$Weight[i] = wts[i]
  }
  return(list("Parm" = Parm, "Q" = Q))
}
#
#
#
#
#
