plotModelDotDensity3D = function(Data, XKernel, YKernel,
                                 EmpiricDataPDE, DotDensity,
                                 Means, Covariances, Weights, MainAxesAngle,
                                 Colors, Cls,
                                 Camera = NULL,
                                 ShowScatter = TRUE,
                                 ShowMarkers = FALSE,
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
  # EmpiricDataPDE[1:n]            Numeric vector with density estimation of
  #                                Data defined for each datapoint within the
  #                                Data.
  # DotDensity[1:n]                Numeric vector with density values for each
  #                                observation from Data.
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
  #if((dim(ContinuousDataPDE)[1] != length(YKernel)) | (dim(ContinuousDataPDE)[2] != length(YKernel))){
  #  message("Parameter ContinuousDataPDE must have dimensions the same size as the length of XKernel. Returning.")
  #  return()
  #}
  if(missing(EmpiricDataPDE)){
    message("Parameter EmpiricDataPDE is missing. Returning.")
    return()
  }else{
    if(!is.vector(EmpiricDataPDE)){
      message("Parameter EmpiricDataPDE is not of type vector. Returning.")
      return()
    }
  }
  if((dim(Data)[1] != length(EmpiricDataPDE))){
    message("Number of rows of parameter Data must match length of vector EmpiricDataPDE. Returning.")
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
  }else{
    if(!is.vector(Cls)){
      message("Parameter Cls is not of type vector. Returning.")
      return()
    }
  }
  if(dim(Data)[1] != length(Cls)){
    message("Number of rows of parameter Data must match length of vector Cls. Returning.")
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
  if(is.null(ShowMarkers)){
    ShowMarkers = FALSE
  }
  if(!is.logical(Debug)){
    message("Parameter Debug is not a logical type. Returning.")
    return()
  }
  if(Debug){
    cat(file = stderr(), "Plot 3D\n")
  }
  EmpiricDataPDE = EmpiricDataPDE/sum(EmpiricDataPDE)
  MinData = min(Data)
  MaxData = max(Data)
  # Compute the density for a density dot plot of the complete GMM model
  plotOut = plotly::plot_ly(source = Source)
  plotOut = plotly::add_trace(p = plotOut,
                              x = Data[,1],
                              y = Data[,2],
                              z = DotDensity,
                              marker = list(color = Colors[Cls], size = 2),
                              type = "scatter3d", mode = "markers")
  #if(ShowMarkers){
  if(FALSE){ # Experimental function deactivated
    for(i in 1:length(Means)){
      plotOut = plotly::add_trace(p = plotOut,
                                  x = Means[[i]][1],
                                  y = Means[[i]][2],
                                  z = MaxDensityPerClass[i],
                                  marker = list(color = "black", size = 10,
                                                line = list(color = 'bisque',
                                                            width = 4)),
                                  type = "scatter3d", mode = "markers")
      Vi   = RetrieveMainAxesInfoFromGMM(Covariances, MainAxesAngle)
      PC1A = Vi$PC1A
      PC1B = Vi$PC1B
      PC2A = Vi$PC2A
      PC2B = Vi$PC2B
      plotOut = plotly::add_trace(p = plotOut,
                                  x = c(Means[[i]][1], Means[[i]][1] + PC1A),
                                  y = c(Means[[i]][2], Means[[i]][2] + PC1B),
                                  z = c(MaxDensityPerClass[i], MaxDensityPerClass[i]),
                                  type = "scatter3d", mode = "lines", name = "",
                                  line = list(width = 8, color = "bisque"))
      plotOut = plotly::add_trace(p = plotOut,
                                  x = Means[[i]][1] + PC1A,
                                  y = Means[[i]][2] + PC1B,
                                  z = MaxDensityPerClass[i],
                                  type = "scatter3d", mode = "markers",
                                  marker = list(size = 6, color = "bisque",
                                                symbol = "152"))
      plotOut = plotly::add_trace(p = plotOut,
                                  x = c(Means[[i]][1], Means[[i]][1] + PC2A),
                                  y = c(Means[[i]][2], Means[[i]][2] + PC2B),
                                  z = c(MaxDensityPerClass[i], MaxDensityPerClass[i]),
                                  type = "scatter3d", mode = "lines", name = "",
                                  line = list(width = 6, color = "bisque"))
      plotOut = plotly::add_trace(p = plotOut,
                                  x = Means[[i]][1] + PC2A,
                                  y = Means[[i]][2] + PC2B,
                                  z = MaxDensityPerClass[i],
                                  type = "scatter3d", mode = "markers",
                                  marker = list(size = 6, color = "bisque",
                                                symbol = "152"))
    }
  }
  # Compute the density for a scatter plot of the Data density estimation
  if(ShowScatter == TRUE){
    plotOut = plotly::add_markers(p = plotOut,
                                  x = Data[,1],
                                  y = Data[,2],
                                  z = EmpiricDataPDE,
                                  marker = list(size = 1.5,
                                                color = "black"))
  }
  plotOut = plotly::event_register(p = plotOut, event = 'plotly_relayout')
  plotOut = plotly::config(p = plotOut, displayModeBar = F)
  plotOut = plotly::hide_colorbar(p = plotOut)
  plotOut = plotly::hide_legend(p = plotOut)
  if(ShowScatter == TRUE){
    Title = "Model PDF (color) and Empirical Dot Density (black)"
  }else{
    Title = "Model PDF"
  }
  plotOut = plotly::layout(p = plotOut,
                           title = Title,
                           scene = list(
                             xaxis = list(title = AxNames[1], range = c(MinData, MaxData)),
                             yaxis = list(title = AxNames[2], range = c(MinData, MaxData)),
                             zaxis = list(title = "Probability"),
                             #zaxis = list(range = c(0, max(EmpiricDataPDE))),
                             aspectratio = list(x = 1, y = 1, z = 1),
                             #scene, Scene))
                             camera = Camera))
  return(plotOut)
}
