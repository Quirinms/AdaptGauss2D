computeRMS = function(Data, Means, Covariances, Weights, EmpiricDataPDE){
  # DESCRIPTION
  # Interface to Expectation Maximization algorithm from Dr. Baggenstoss.
  #
  # INPUT
  # Data[1:n, 1:2]                 Numeric matrix with n observations and 2 features.
  # Means                          List with l [1:2] numerical vector defining
  #                                the means of the l GMM components.
  # Covariances                    List with l [1:2, 1:2] numerical matrices
  #                                defining the covariance matrices of the l GMM
  #                                components.
  # Weights[1:l]                   Numerical vector with weights for each GMM
  #                                component.
  # EmpiricDataPDE[1:n]            Numeric vector with density estimation of
  #                                Data defined for each datapoint within the
  #                                Data.  #
  # OUTPUT
  # Res    List of GMM settings computed by an EM algorithm of choice.
  #
  #
  # Author: QMS 15.12.2021

  # Compute the density for a density dot plot
  GMMDensity = 0
  MaxDensityPerClass = c()
  for(i in 1:length(Means)){
    TmpDensity         = mixtools::dmvnorm(y     = Data,               # Marginal density estimation
                                           mu    = Means[[i]],
                                           sigma = Covariances[[i]])
    # Ensure probability density as part of mixture (respect weights of each component)
    TmpDensity         = Weights[i] * TmpDensity/sum(TmpDensity)
    MaxDensityPerClass = c(MaxDensityPerClass, max(TmpDensity))
    GMMDensity         = GMMDensity + TmpDensity
  }
  RMS = sqrt(sum((GMMDensity-EmpiricDataPDE)**2))*100

  return(round(RMS, 5))
}
