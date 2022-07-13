plotScatter = function(Data, CurrGauss,
                       Means, Covariances, MainAxesAngle,
                       Colors, Cls,
                       Shapes, ShapeText, AxNames = c("X1", "X2"),
                       ShowAxis = FALSE, ShowEllipsoids = FALSE,
                       ShowGaussNr = FALSE, Source = "D"){
  # DESCRIPTION
  #
  #
  # INPUT
  # Data                  [1:n, 1:2] Numeric matrix with n observations and 2 features.
  # Means                 List with l [1:2] numerical vector defining the means of
  #                       the l GMM components.
  # Covariances           List with l [1:2, 1:2] numerical matrices defining the
  #                       covariance matrices of the l GMM
  #                       components.
  # MainAxesAngle[1:4]    Numeric vector with 1st and 2nd main axes of a 2D
  #                       ellipsoid and the respective angles
  #                       measured to the first unit vector c(0,1).
  # Colors                [1:n] Numerical vector of size n containing the colors of each observation.
  # Cls                   [1:n] Numerical vector of size n containing the classes of each observation.
  # Shapes                List of List with 4 attributes (type, fillcolor, opacity, path)
  #                       for a shape for plotting. Here it is used for plotting an ellipsoid.
  # AxNames               Character vector with names for each dimension ax of 2D
  #                       plot.
  # ShowAxis              Boolean (Default=TRUE). T: Show main axes of the GMM
  #                       ellipsoid/circle.
  # Source                Character indicating plot source. Important attribute for plotly in
  #                       shiny in order to keep control of specific panels.
  #
  # OUTPUT
  # plotOut    Plotly object containing plot for direct visualization.
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
  }
  if(dim(Data)[1] != length(Cls)){
    message("Number of rows of parameter Data must match length of vector Cls. Returning.")
    return()
  }
  if(length(Colors) < length(unique(Cls))){
    message("Length of parameter Colors must be greater than or equal to the number of unique entries in Cls. Returning.")
    return()
  }
  if(!is.list(Shapes)){
    message("Parameter Shapes must be of type list. Returning.")
    return()
  }
  if(is.null(ShowAxis)){
    ShowAxis = TRUE
  }
  if(!is.character(Source)){
    message("Parameter Source is not a character type. Returning.")
    return()
  }
  plotOut = plotly::plot_ly(source = Source)
  plotOut = plotly::add_markers(p = plotOut,
                                x = Data[,1],
                                y = Data[,2],
                                marker = list(color = Colors[Cls]), type = "scatter")
  if(ShowGaussNr){
    plotOut = plotly::add_text(p = plotOut,
                               x = ShapeText[,1],
                               y = ShapeText[,2],
                               text  = ShapeText[,3],
                               textfont = list(color = "white", size = 40))
  }
  if(ShowAxis){
    for(i in 1:length(Means)){
      plotOut = plotly::add_markers(p = plotOut,
                                    x = Means[[i]][1],
                                    y = Means[[i]][2],
                                    marker = list(color = "white", size = 7,
                                                  line = list(color = "black", width = 2)),
                                    type = "scatter")
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
                                          arrowcolor = "bisque", arrowhead = 0.7, arrowsize = 2)
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
  plotOut = plotly::layout(p      = plotOut,
                           title  = "2D Scatter Plot of Dataset",
                           shapes = Shapes,
                           xaxis  = list(title = AxNames[1], fixedrange = T, scaleanchor="y", scaleratio=1),
                           yaxis  = list(title = AxNames[2], fixedrange = T),
                           plot_bgcolor = "rgb(254, 254, 254)",              # plot_bgcolor = "rgb(254, 247, 234)",
                           paper_bgcolor = "rgb(254, 254, 254)")             # paper_bgcolor = "rgb(254, 247, 234)"
  plotOut = plotly::hide_colorbar(p = plotOut)
  plotOut = plotly::hide_legend(p = plotOut)
  plotOut = plotly::config(p = plotOut, displayModeBar=F, editable=T)
  return(plotOut)
}
