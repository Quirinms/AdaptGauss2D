computeGMMClassification = function(Data, Means, Covariances, Weights){
  #
  # DESCRIPTION
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
  # Classification[1:n]          Numeric vector with cls for each data point.
  # CurrDataDensity[1:n, 1:l]    Numeric matrix with probability density for
  #                              all n data points and l classes
  #
  # Author QS 2021

  if(missing(Data)){
    message("Parameter Data is missing. Returning.")
    return()
  }else{
    if(!is.matrix(Data)){
      message("Parameter Data is not of type matrix. Returning.")
      return()
    }else if(dim(Data)[2] != 2){
      message("Parameter Data does not have exactly two feature columns. Returning.")
      return()
    }
  }

  if(!is.null(Means)){
    if(!is.list(Means)){
      message("Parameter Means is not of type list. Returning.")
      return()
    }else{
      for(i in 1:length(Means)){
        if(!is.vector(Means[[i]])){
          message("Parameter Means can only contain vectors. Returning.")
          return()
        }else if(length(Means[[i]]) != 2){
          message("Parameter Means can only contain vectors of dimension 2. Returning.")
          return()
        }
      }
    }
  }

  if(!is.null(Covariances)){
    if(!is.list(Covariances)){
      message("Parameter Cov is not of type list. Returning.")
      return()
    }else{
      for(i in 1:length(Covariances)){
        if(!is.matrix(Covariances[[i]])){
          message("Parameter Cov can only contain matrices. Returning.")
          return()
        }else if((dim(Covariances[[i]])[1] != 2) | (dim(Covariances[[i]])[2] != 2)){
          message("Parameter Cov can only contain matrices of dimension 2x2. Returning.")
          return()
        }
      }
    }
  }

  if(!is.null(Weights)){
    if(!is.vector(Weights)){
      message("Parameter Weights is not of type vector. Returning.")
      return()
    }else if(!is.numeric(Weights[i])){
      message("Parameter Weights can only contain numerics. Returning.")
      return()
    }
  }

  gDen = sapply(1:length(Means), function(i){       # density for each point and gaussian
    mixtools::dmvnorm(y = Data, mu = Means[[i]], sigma = Covariances[[i]]) * Weights[i]
  })
  Classification  = apply(gDen, 1, which.max)
  return(list("Classification"  = Classification,
              "CurrDataDensity" = gDen))
}
#
#
#
#
#
