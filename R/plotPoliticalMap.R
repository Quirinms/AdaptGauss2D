plotPoliticalMap = function(Data, CurrGauss,
                            Means, Covariances, Weights, MainAxesAngle,
                            Colors, Cls,
                            Shapes, ShapeText, AxNames = c("X1", "X2"),
                            ShowAxis = FALSE, ShowEllipsoids = TRUE,
                            ShowScatter = FALSE, Source = "D"){
  # DESCRIPTION
  # Draws a political map of the classification of a Gaussian Mixture Model.
  #
  # INPUT
  # Data[1:n, 1:2]    Numeric matrix with n observations and 2 features.
  # Means             List with l [1:2] numerical vector defining the means of
  #                   the l GMM components.
  # Covariances       List with l [1:2, 1:2] numerical matrices defining the
  #                   covariance matrices of the l GMM components.
  # Weights[1:l]      Numerical vector with weights for each GMM component.
  # MainAxesAngle     List of numeric vectors with 1st and 2nd main axes
  #                   of a 2D ellipsoid and the respective angles
  #                   measured to the first unit vector c(0,1).
  # Colors[1:n]       Numerical vector of size n containing the colors of each
  #                   observation.
  # Cls[1:n]          Numerical vector of size n containing the classes of each
  #                   observation.
  # Shapes            List of List with 4 attributes (type, fillcolor, opacity, path)
  #                   for a shape for plotting. Here it is used for plotting an ellipsoid.
  # ShapeText         [1:l, 1:3] Numeric matrix with l means and two entries for the
  #                   two-dimensional coordinates and one entry for the number of the Gaussian component.
  # AxNames           Character vector with names for each dimension ax of 2D
  #                   plot.
  # ShowAxis          Boolean indicating if main axis of models components are
  #                   shown (TRUE) or not (FALSE).
  # ShowEllipsoids    Boolean indicating if ellipsoids of models components are
  #                   shown (TRUE) or not (FALSE).
  # ShowScatter       Boolean indicating if 3D scatter points are shown in plot
  #                   (TRUE) or not (FALSE).
  # Source             Character indicating plot source (Default = "D"). Important
  #                    attribute for plotly in shiny in order to keep control of specific panels.
  #
  # OUTPUT
  # plotOut    Plotly object containing plot for direct visualization.
  #
  # Author: QMS 03.01.2022
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
      message("Parameter CovMatrices is not of type list. Returning.")
      return()
    }else{
      for(i in 1:length(Covariances)){
        if(!is.matrix(Covariances[[i]])){
          message("Parameter CovMatrices can only contain matrices. Returning.")
          return()
        }else if((dim(Covariances[[i]])[1] != 2) | (dim(Covariances[[i]])[2] != 2)){
          message("Parameter CovMatrices can only contain matrices of dimension 2x2. Returning.")
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
  if(missing(MainAxesAngle)){
    message("Parameter MainAxesAngle is missing. Returning.")
    return()
  }else{
    if(!is.list(MainAxesAngle)){
      message("Parameter MainAxesAngle is not of type list. Returning.")
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
  if(!is.logical(ShowScatter)){
    message("Parameter Show3DPoints is not a logical type. Returning.")
    return()
  }
  df = as.data.frame(cbind(Data, as.vector(Cls)))
  df[,3] = as.factor(df[,3])
  grid = as.matrix(base::expand.grid(seq(min(df[, 1]), max(df[, 1]), length.out=100),
                           seq(min(df[, 2]), max(df[, 2]), length.out=100)))
  gDen = sapply(1:length(Means), function(i){       # density for each point and gaussian
    mixtools::dmvnorm(y = grid, mu = Means[[i]], sigma = Covariances[[i]]) * Weights[i]
  })
  matchingGauss = apply(gDen, 1, which.max)
  plotOut = plotly::plot_ly(data = data.frame(grid),
                            x = grid[,1], y = grid[,2],
                            mode = "markers",
                            type = "scatter",
                            color = matchingGauss,
                            colors = Colors[1:length(unique(matchingGauss))],
                            source = Source)
  if(ShowScatter){
    plotOut = plotly::add_markers(p = plotOut,
                                  x = Data[,1],
                                  y = Data[,2],
                                  color = Colors[Cls[]],
                                  marker = list(size = 3, color = "black"))#, colors = Colors[1:length(unique(Cls))])
  }
  if(ShowAxis){
    for(i in 1:length(Means)){
      plotOut = plotly::add_markers(p = plotOut,
                                    x = Means[[i]][1],
                                    y = Means[[i]][2],
                                    marker = list(color = "bisque"), type = "scatter")
      if(all(Covariances[[i]] != diag(c(1,1)))){
        MySVD = svd(Covariances[[i]])                                   # Compute singular value decomposition for Princ. Component Axes
        SD1 = MySVD$d[1]*MySVD$u[,1]                                 # Extract 1st PCA component vector
        SD2 = MySVD$d[2]*MySVD$u[,2]                                 # Extract 2nd PCA component vector
        NormSD1   = norm(SD1, type = "2")
        TopCircle1    = acos(sum(SD1 * c(0,1))/NormSD1)*(180/pi)  # See if 1st PCA is on the upper part of the cartesian coord. sys.
        BottomCircle1 = acos(sum(SD1 * c(0,-1))/NormSD1)*(180/pi) # See if 1st PCA is on the lower part of the cartesian coord. sys.
        Angle1 = acos(sum(SD1 * c(1,0))/NormSD1)*(180/pi)
        if(BottomCircle1<TopCircle1){
          Angle1 = 360 - Angle1                                   # This would be the angle for the lower part
        }
        if(round(abs(Angle1-MainAxesAngle[[i]][3])) > 5 & MainAxesAngle[[i]][3] != 360){
          SD1 = -SD1
        }
        NormSD2   = norm(SD2, type = "2")
        TopCircle2    = acos(sum(SD2 * c(0,1))/NormSD2)*(180/pi) # See if 1st PCA is on the upper part of the cartesian coord. sys.
        BottomCircle2 = acos(sum(SD2 * c(0,-1))/NormSD2)*(180/pi)
        Angle2 = acos(sum(SD2 * c(1,0))/NormSD2)*(180/pi)
        if(BottomCircle2<TopCircle2){
          Angle2 = 360 - Angle2                                   # This would be the angle for the lower part
        }
        if(abs(Angle2-((MainAxesAngle[[i]][3]+90)%%360)) > 5){
          SD2 = -SD2
        }
        PC1A = SD1[1]; PC1B = SD1[2]; PC2A = SD2[1]; PC2B = SD2[2] # Eigenvector components
        plotOut = plotly::add_annotations(p = plotOut,
                                          standoff=0,
                                          x = Means[[i]][1] + PC1A, y = Means[[i]][2] + PC1B,
                                          ax = Means[[i]][1], ay = Means[[i]][2],
                                          xref = "x", yref = "y",
                                          axref = "x", ayref = "y",
                                          text = "", showarrow = TRUE,
                                          arrowcolor="bisque", arrowhead = 0.7, arrowsize = 2)
        #plotOut = plotly::add_annotations(p = plotOut,
        #                                  x = Means[[i]][1] + PC2A, y = Means[[i]][2] + PC2B,
        #                                  ax = Means[[i]][1], ay = Means[[i]][2],
        #                                  xref = "x", yref = "y",
        #                                  axref = "x", ayref = "y",
        #                                  text = "", showarrow = TRUE,
        #                                  arrowcolor="bisque", arrowhead = 0.7, arrowsize = 1)
      }
    }
  }
  #minData = round(min(min(Data[,1]), min(Data[,2])), 2)
  #minData = minData + sign(minData) * 0.1 * abs(minData)
  #maxData = round(max(max(Data[,1]), max(Data[,2])), 2)
  #maxData = maxData + sign(maxData) * 0.1 * abs(maxData)
  Xaxis <- list(title = AxNames[1],
                fixedrange = T, scaleanchor="y", scaleratio=1,
                zeroline = FALSE,
                showline = FALSE,
                showticklabels = TRUE,
                showgrid = FALSE)
  Yaxis <- list(title = AxNames[2],
                fixedrange = T,
                zeroline = FALSE,
                showline = FALSE,
                showticklabels = TRUE,
                showgrid = FALSE)
  if(length(Shapes) == length(Means)){
    for(i in 1:length(Means)){
      if(i != CurrGauss){
        Shapes[[i]]$fillcolor = "black"
      }else{
        Shapes[[i]]$opacity = 0.7
      }
    }
  }
  if(ShowEllipsoids != TRUE){
    Shapes = NULL
  }
  plotOut = plotly::layout(p = plotOut,
                           shapes = Shapes,
                           title = "Political Map",
                           xaxis = Xaxis,
                           yaxis = Yaxis)
  plotOut = plotly::hide_colorbar(p = plotOut)
  plotOut = plotly::hide_legend(p = plotOut)
  plotOut = plotly::config(p = plotOut, displayModeBar=F, editable=T)
  return(plotOut)
}
