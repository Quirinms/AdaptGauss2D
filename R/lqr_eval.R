lqr_eval = function(Data, ModeMean, CholeskyCovar){
  # V = lqr_eval(Data, ModeMean, CholeskyCovar)
  # V$logPDF
  #
  # calculates the Gaussian logPDF given Mean and Covariances
  # Input data is dimensioned (DIM-by-N)
  #
  # INPUT
  # Data[1:N,1:d]             Numerical matrix, N = number of cases, d = dimension
  # ModeMean[1:d]             Numerical vector with one d-dim mode in GMM
  # CholeskyCovar[1:d,1:d]    Numerical matrix with covariance matrix of d-dim mode
  #                           (dxd matrix where d denotes the feature's dimension).
  #
  # OUTPUT
  # logPDF[1:N]               Numerical vector with log Gaussian PDF.
  #
  # in /dbt/EMforGauss/
  # Author QS 2021
  if(is.null(Data)){
    warning("Input parameter Data is not given. Returning.")
    return()
  }
  if(is.null(ModeMean)){
    warning("Input parameter ModeMean is not given. Returning.")
    return()
  }
  if(is.null(CholeskyCovar)){
    warning("Input parameter CholeskyCovar is not given. Returning.")
    return()
  }
  if(!is.matrix(ModeMean)){
    ModeMean = as.matrix(ModeMean)
  }
  if(dim(ModeMean)[1] != 1){ # Assert vertical vector
    ModeMean = t(ModeMean)
  }
  N   = dim(Data)[1]
  DIM = dim(Data)[2]
  DIM0= length(ModeMean) # Ueberpruefung ob das Modell auf die Daten passt
  if(DIM0 != DIM){
    warning("lqr_eval: data and given Mean have different dimension")
  }
  if(dim(CholeskyCovar)[1] != dim(CholeskyCovar)[2]){
    warning("lqr_eval: CholeskyCovar not square matrix.")
  }
  if(dim(CholeskyCovar)[1] != DIM){
    warning("lqr_eval: dimensions of data and CholeskyCovar does not match.")
  }
  AbsDiag = abs(diag(CholeskyCovar))
  if(sum(AbsDiag) == 0){ # abfangen log(0)
    l = matrix(0,N,1)
  }else{
    logdet = 2 * sum(log(AbsDiag))
    tmpidx = Data - t(matrix(ModeMean, length(ModeMean), N))
    tmpidx = 0.5 * colSums((pracma::mldivide(t(CholeskyCovar), t(tmpidx)))**2)
    l = -tmpidx - (DIM/2)*log(2*pi) - 0.5*logdet
  }
  return(list("LogPDF"=l))
}
#
#
#
#
#
