computeGridDensity = function(XKernel, YKernel,
                              Means, Covariances, Weights){
  # DESCRIPTION
  # Interface to Expectation Maximization algorithm from Dr. Baggenstoss.
  #
  # INPUT
  # XKernel[1:x]      Numeric vector defining domain of x axis.
  # YKernel[1:x]      Numeric vector defining domain of y axis.
  # Means             List with l [1:2] numerical vector defining the means of
  #                   the l GMM components.
  # Covariances       List with l [1:2, 1:2] numerical matrices defining the
  #                   covariance matrices of the l GMM
  #                   components.
  # Weights[1:l]      Numerical vector with weights for each GMM component.
  #
  # OUTPUT
  # GridDensity[1:x, 1:x]    Numeric matrix PDF for a GMM with one or more
  #                          Gaussians on a Grid defined by the vectors XKernel
  #                          and YKernel.
  #
  # Author: QMS 13.04.2021
  DomainGrid = as.matrix(expand.grid(XKernel, YKernel))
  GridDensity = 0
  for(i in 1:length(Means)){
    TmpDensity = mvtnorm::dmvnorm(x     = DomainGrid,        # Marginal density estimation
                                  mean  = Means[[i]],
                                  sigma = Covariances[[i]])
    # Weight each componente with its respective weight and ensure total sum = 1
    TmpDensity = Weights[i] * TmpDensity/sum(TmpDensity)
    GridDensity = GridDensity + TmpDensity
  }
  GridDensity = matrix(data = GridDensity, nrow = length(XKernel), ncol = length(YKernel), byrow = T)
  return(GridDensity)
}
