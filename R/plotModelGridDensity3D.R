plotModelGridDensity3D = function(Data,
                                  XKernel, YKernel,
                                  ContinuousDataPDE,
                                  GridDensity,
                                  Means, Covariances, Weights,
                                  MainAxesAngle,
                                  Colors, Cls,
                                  Camera = NULL,
                                  AxNames = c("X", "Y"),
                                  Source = "F1", Debug = FALSE){
  # DESCRIPTION1
  # Compares the continuous density plot of a dataset with the density
  # distribution defined by one or more gaussians (GMM; Mean, Cov, Weights).
  # The datapoints from the dataset have the density assigned by the GMM
  # component with the the highest probability (class assignment must be given
  # and is not computed here).
  #
  # INPUT
  # Data[1:n, 1:2]                 Numeric matrix with n observations and 2
  #                                features.
  # XKernel[1:x]                   Numeric vector defining domain of x axis.
  # YKernel[1:x]                   Numeric vector defining domain of y axis.
  # GridDensity[1:x, 1:x]          Numeric matrix PDF for a GMM with one or more
  #                                Gaussians on a Grid defined by the vectors XKernel
  #                                and YKernel.
  # Means                          List with l [1:2] numerical vector defining
  #                                the means of the l GMM components.
  # Covariances                    List with l [1:2, 1:2] numerical matrices
  #                                defining the covariance matrices of the l GMM
  #                                components.
  # Weights[1:l]                   Numerical vector with weights for each GMM
  #                                component.
  # MainAxesAngle[1:4]             Numeric vector with 1st and 2nd main axes
  #                                of a 2D ellipsoid and the respective angles
  #                                measured to the first unit vector c(0,1).
  # Colors[1:n]                    Numerical vector of size n containing the
  #                                colors of each observation.
  # Cls[1:n]                       Numerical vector of size n containing the
  #                                classes of each observation.
  # Camera                         List of attributes concerning the camera
  #                                angle for plotly visualizations.
  # ShowScatter                    Boolean (Default=TRUE). T: Show Data
  #                                Scatter with empirical estimated density.
  #                                F: Do not show nothing.
  # ShowMarkers                    Boolean (Default=TRUE). T: Show main axes
  #                                of the GMM ellipsoid/circle.
  # AxNames                        Character vector with names for each
  #                                dimension ax of 2D plot.
  # Source                         (Default = "F1"). Character indicating plot
  #                                source. Important attribute for plotly in
  #                                shiny in order to keep control of specific
  #                                panels.
  # Debug                          Boolean (Default=FALSE). T: Show developer
  #                                information and warnings in terminal.
  #                                F: Show nothing.
  #
  # OUTPUT
  # plotOut    Plotly object containing plot for direct visualization.
  #
  # DETAILS
  # Uses a density estimation (discretized computation on grid) to display
  # a continuous density plot on domain defined by XKernel and YKernel of the
  # whole empircal data distribution.
  # The density of the GMM components is computed via an multivariate gaussian
  # formula.
  # The datapoints take the density of the GMM components with highest
  # probability to take the specific datapoint. This computation must be done
  # outside of this routine (give the Cls/Class/Labels as parameter).
  # Parameter Camera is important for shiny to let the view fixed during
  # continuous updates of the interactive visualization (Shiny App). The Source
  # is important to separate between multiple visualizations in the Shiny App.
  #
  # Author: QMS 15.12.2021
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
  if(missing(XKernel)){
    message("Parameter XKernel is missing. Returning.")
    return()
  }else{
    if(!is.vector(XKernel)){
      message("Parameter XKernel is not of type vector. Returning.")
      return()
    }
  }
  if(missing(YKernel)){
    message("Parameter YKernel is missing. Returning.")
    return()
  }else{
    if(!is.vector(YKernel)){
      message("Parameter YKernel is not of type vector. Returning.")
      return()
    }
  }
  if(length(XKernel) != length(YKernel)){
    message("Parameter XKernel and YKernel must be vectors of same length. Returning.")
    return()
  }
  if((dim(ContinuousDataPDE)[1] != length(YKernel)) | (dim(ContinuousDataPDE)[2] != length(YKernel))){
    message("Parameter ContinuousDataPDE must have dimensions the same size as the length of XKernel. Returning.")
    return()
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
  if(missing(Colors)){
    message("Parameter Colors is missing. Returning.")
    return()
  }else{
    if(!is.vector(Colors)){
      message("Parameter Colors is not of type vector. Returning.")
      return()
    }
  }
  if(missing(Cls)){
    message("Parameter Cls is missing. Returning.")
    return()
  }
  if(length(Colors) < length(unique(Cls))){
    message("Length of parameter Colors must be greater than or equal to the number of unique entries in Cls. Returning.")
    return()
  }
  if(!is.character(Source)){
    message("Parameter Source is not a character type. Returning.")
    return()
  }
  if(!is.logical(Debug)){
    message("Parameter Debug is not a logical type. Returning.")
    return()
  }
  if(Debug){
    cat(file = stderr(), "Plot 3D\n")
  }
  MinData = min(Data)
  MaxData = max(Data)
  #DomainGrid = as.matrix(expand.grid(XKernel, YKernel))
  # Compute the density for a density dot plot
  #GMMDensity = 0
  #MaxDensityPerClass = c()
  #for(i in 1:length(Means)){
  #  TmpDensity         = mvtnorm::dmvnorm(y = DomainGrid,               # Marginal density estimation
  #                                 mu = Means[[i]],
  #                                 sigma = Covariances[[i]])
  #  # Ensure probability density as part of mixture (respect weights of each component)
  #  TmpDensity         = Weights[i] * TmpDensity/sum(TmpDensity)
  #  MaxDensityPerClass = c(MaxDensityPerClass, max(TmpDensity))
  #  GMMDensity         = GMMDensity + TmpDensity
  #}
  #GMMDensity = GMMDensity/sum(GMMDensity)
  #GMMDensity = t(matrix(data = GMMDensity, nrow = length(XKernel), ncol = length(YKernel)))
  colfunc <- colorRampPalette(c("white","black"))
  MyColorGradient = colfunc(10)
  if(requireNamespace("colorRamps"))
	MyColorGradient = colorRamps::matlab.like(10)
  else
	MyColorGradient = colfunc(10)
  #MyColorGradient = c(rep("#FFFFFF", 2), MyColorGradient)
  plotOut = plotly::plot_ly(x = XKernel,
                            y = YKernel,
                            z = ContinuousDataPDE,
                            type = "surface",
                            colors = MyColorGradient,
                            alpha = 1)
  plotOut = plotly::add_surface(p          = plotOut,
                                z          = GridDensity,
                                colorscale = list(c(0,1), c("rgb(107,184,255)","rgb(0,90,124)")),
                                opacity    = 0.6)
  plotOut = plotly::hide_colorbar(p = plotOut)
  plotOut = plotly::layout(p          = plotOut,
                           title      = "Model (blue) vs. Emp. Grid Dens.",
                           showlegend = FALSE,
                           scene      = list(
                             xaxis = list(title = AxNames[1], range = c(MinData, MaxData)),
                             yaxis = list(title = AxNames[2], range = c(MinData, MaxData)),
                             zaxis = list(title = "Probability"),
                             #zaxis = list(range = c(0, max(EmpiricDataPDE))),
                             aspectratio = list(x = 1, y = 1, z = 1),
                             #scene, Scene))
                             camera = Camera))
  return(plotOut)
}
