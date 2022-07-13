computeDotDensity = function(Data, Means, Covariances, Weights){
  # DESCRIPTION
  #
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
  #
  # OUTPUT
  # DotDensity[1:n]    Numeric vector with density values for each observation
  #                    from Data.
  #
  # Author: QMS 15.12.2021
  # Compute the density for a density dot plot
  DotDensity = 0
  for(i in 1:length(Means)){
    TmpDensity = mixtools::dmvnorm(y     = Data,               # Marginal density estimation
                                   mu	 = Means[[i]],
                                   sigma = Covariances[[i]])
    # Ensure probability density as part of mixture (respect weights of each component)
    TmpDensity = Weights[i] * TmpDensity/sum(TmpDensity)
    DotDensity = DotDensity + TmpDensity
  }
  DotDensity = DotDensity/sum(DotDensity)
  return(DotDensity)
}
