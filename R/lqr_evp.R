lqr_evp = function(Parm, Data, Flag){
  # V = lqr_evp(Parm, Data, Flag)
  # V$lg
  # V$ModePDFs
  # V$
  # PDF
  # DESCRIPTION
  #
  # Returns either total log-likelihood of Gaussian mixtures or separate log
  # likelihoods for all modes.
  #
  #
  # Input
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
  # Flag             Boolean Flag
  #                          Flag=1,
  #                          if(  Flag==1, returns log mode PDF in columns, that
  #                          is the output of lqr_eval only, not including
  #                          mixing weights for each column.
  #                          Flag=0
  #                          Computes the total log-PDF for normalized data in a
  #                          single output column (i.e. the log of the weighted
  #                          sum of the rows of the Flag=1 output) plus the
  #                          Jacobian (it outputs the  PDF value with respect to
  #                          the raw unnormalized data by taking into account
  #                          the jacobian of the normalization operation).
  #                          In this case outputs ModePDFs and PDF are derived.
  #
  # OUTPUT
  # lg[1:n]              Numerical vector with LogPDF
  # ModePDFs[1:n,1:L]    Numerical matrix with PDFs for all single modes
  # PDF[1:n]             Numerical vector with combined Gaussian as weighthed
  #                      sum of PDF
  #
  #  derived from : Dr. Paul M. Baggenstoss
  #  optimized by ALU
  nmode = length(Parm$Modes$Weight) # the number of modes
  wts   = Parm$Modes$Weight # the weights of the modes, only used for (Flag==0)
  DIM0  = length(Parm$Features$Names) # dimensionality of the GMMs
  N     = dim(Data)[1]
  DIM   = dim(Data)[2]
  if(DIM != DIM0) {
    warning('lqr_evp: data and parms have different dimensions')
  }
  lg = matrix(0, N,nmode)
  for(k in 1:nmode){  # ausrechnen der LogPDF fuer den k-ten Modus
    idxFrom        = DIM*(k-1)+1
    idxTo          = DIM*k
    cholesky_covar = Parm$Modes$Covariance[idxFrom:idxTo,]
    V              = lqr_eval(Data, Parm$Modes$Mean[k,], cholesky_covar)
    lg[,k]         = V$LogPDF
  }
  # damit is Flag ==1 erledigt:
  # lg== log mode PDF in columns, that is the  output of lqr_eval(...)
  # not including mixing weights for each column.
  if(Flag==0){# lg ==  total log-PDF for normalized data in a single  output column
  #   i.e. the log of the weighted sum of the rows of lg
  #   plus the Jacobian (it outputs the
  #   PDF value with respect to the raw unnormalized data by
  #   taking into account the jacobian of the normalization operation).
    if(nmode > 1){
      # Max row of matrix
      mx = apply(lg, 2, max)    # mx == max of each  mode PDF for each data point
      Allmax = max(mx)
      for(i in 1:nmode){  # subtract that from each column        => max(LogPDF(:,i)) == 0
        lg[,i] = lg[,i] - mx[i]
      }
      ModePDFs = exp(lg)  # compute the PDFs for all single modes => max(ModePDFs(:,i)== 1
      for(i in 1:nmode){
        ModePDFs[,i] = ModePDFs[,i] * wts[i]/Allmax  # weighting of the ModePDFS
      }
      PDF = rowSums(ModePDFs)#  Compute combined Gaussian as weighthed sum of PDF
      lg  = log(PDF)          # lg = log(PDF) + mx == add back in the scaling ???
    }
  }else{
    ModePDFs = NULL     # no output
    PDF      = NULL
  }
  return(list("LogPDF"=lg, "ModePDFs"=ModePDFs, "PDF"=PDF))
}
