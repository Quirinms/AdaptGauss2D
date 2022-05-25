plotDensity = function(Data, Cls, CurrGauss, Colors,
                       GridDensity,
                       Means, Covariances, Weights, MainAxesAngle,
                       XKernel, YKernel,
                       Shapes, ShapeText, AxNames = c("X1", "X2"),
                       ShowAxis = FALSE, ShowEllipsoids = TRUE,
                       ShowScatter = FALSE,
                       ShowGaussNr = FALSE, Source = "D"){
  # DESCRIPTION
  #
  #
  # INPUT
  # Data[1:n, 1:2]     Numeric matrix with n observations and 2 features.
  # Cls[1:n]           Numerical vector of size n containing the classes of each
  #                    observation.
  # Means              List with l [1:2] numerical vector defining the means of
  #                    the l GMM components.
  # Covariances        List with l [1:2, 1:2] numerical matrices defining the
  #                    covariance matrices of the l GMM
  #                    components.
  # Weights[1:l]       Numerical vector with weights for each GMM component.
  # MainAxesAngle[1:4] Numeric vector with 1st and 2nd main axes of a 2D
  #                    ellipsoid and the respective angles
  #                    measured to the first unit vector c(0,1).
  # XKernel[1:x]       Numeric vector defining domain of x axis.
  # YKernel[1:x]       Numeric vector defining domain of x axis.
  # Shapes             List of List with 4 attributes (type, fillcolor, opacity, path)
  #                    for a shape for plotting. Here it is used for plotting an ellipsoid.
  # ShapeText          [1:l, 1:3] Numeric matrix with l means and two entries for the
  #                    two-dimensional coordinates and one entry for the number of the Gaussian component.
  # AxNames            Character vector with names for each dimension ax of 2D
  #                    plot.
  # ShowAxis           Boolean (Default=TRUE). T: Show main axes of the GMM
  #                    ellipsoid/circle.
  # ShowGaussNr        Boolean determining if the shown Gaussian components are
  #                    enumerated (TRUE) or not (FALSE).
  # Source             Character indicating plot source (Default = "D"). Important
  #                    attribute for plotly in shiny in order to keep control of specific panels.
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

  if(missing(Cls)){
    message("Parameter Cls is missing. Returning.")
    return()
  }

  if(!is.null(Means)){
    if(!is.list(Means)){
      message("Parameter Means is not of type list. Returning.")
      return()
    }else{
      if(length(Means) == 0){
        message("Parameter Means is of length 0. Returning.")
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

  if(!is.list(Shapes)){
    message("Parameter Shapes must be of type list. Returning.")
    return()
  }

  if(!is.null(ShapeText)){
    if(!is.matrix(ShapeText)){
      message("Parameter ShapeText is not a logical type. Returning.")
      return()
    }
  }

  if(!is.logical(ShowGaussNr)){
    message("Parameter ShowGaussNr is not a logical type. Returning.")
    return()
  }

  if(!is.character(Source)){
    message("Parameter Source is not a character type. Returning.")
    return()
  }

  #MyPalette = colorRampPalette(colors = c("white", "green"))(2)

  # Compute the density for a continuous plot
  #TestGrid = as.matrix(expand.grid(XKernel, YKernel))
  #GMMDensity = 0
  #for(i in 1:length(unique(Cls))){
  #  #TmpIdx = which(Cls == i)                              # Process each class
  #  TmpDensity = mvtnorm::dmvnorm(x     = TestGrid,       # Marginal density estimation
  #                                 mu    = Means[[i]],
  #                                 sigma = Covariances[[i]])
  #  TmpDensity = Weights[i] * TmpDensity/sum(TmpDensity)   # Ensure probability density as part of mixture (respect weights of each component)
  #  GMMDensity = GMMDensity + TmpDensity
  #}
  #GridDensity = matrix(GMMDensity, nrow = length(XKernel), ncol = length(XKernel), byrow = T)

  # Colorgradient from white to blue
  #colfunc = colorRampPalette(c("white", "blue"))
  #MyColorGradient = c("#FFFFFF", colfunc(10))
  #MyColorGradient = c("#FFFFFF", "#FFFFFF", rainbow(9))
  # Colorgradient rainbow colors
  #MyColorGradient = colorRamps::matlab.like(5) # do not choose to many colors, so that 0 density values are displayed white
  #MyColorGradient = c("#FFFFFF", "#FFFFFF", MyColorGradient) # White out the first values so that the background appears white
  #c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF")
  # Compute the density for a continuous plot of the complete GMM model
  plotOut = plotly::plot_ly(x = XKernel,
                            y = YKernel,
                            type = "contour",
                            colors = PmatrixColormap,
                            source = Source)
  plotOut = plotly::add_contour(p = plotOut,
                                z = GridDensity)

  # Use the Data to estimate a Data Density dot plot
  # Color scales: "Blues", "bilbao"
  #plotOut = plotly::plot_ly(source = Source)
  #plotOut = plotly::add_trace(p = plotOut, x=Data[,1], y  =  Data[,2],
  #                            type = "histogram2dcontour", reversescale = F, color = "bilbao")#,
  #                            #contours = list(start = 0, end = 100, size = 5),
  #                            #colorscale = MyPalette)
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

  if(length(Shapes) == length(Means)){
    for(i in 1:length(Means)){
      if(i != CurrGauss){
        Shapes[[i]]$fillcolor = "black"
      }else{
        Shapes[[i]]$opacity = 0.7
      }
    }
  }

  if(ShowScatter){
    plotOut = plotly::add_markers(p = plotOut,
                                  x = Data[,1],
                                  y = Data[,2],
                                  color = Colors[Cls[]],
                                  marker = list(size = 3, color = "black"))#, colors = Colors[1:length(unique(Cls))])
  }

  if(ShowEllipsoids != TRUE){
    Shapes = NULL
  }

  plotOut = plotly::layout(p      = plotOut,
                           title  = "2D Model Density",
                           shapes = Shapes,
                           xaxis = list(title = AxNames[1],
                                        range = c(min(Data[,1], max(Data[,1]))),
                                        fixedrange = F,
                                        scaleanchor="y",
                                        scaleratio=1,
                                        zeroline = FALSE,
                                        showline = FALSE,
                                        showticklabels = T,
                                        showgrid = FALSE),
                           yaxis = list(title = AxNames[2],
                                        range = c(min(Data[,2], max(Data[,2]))),
                                        fixedrange = F,
                                        zeroline = FALSE,
                                        showline = FALSE,
                                        showticklabels = T,
                                        showgrid = FALSE),
                           plot_bgcolor = "rgb(254, 254, 254)",
                           paper_bgcolor = "rgb(254, 254, 254)")
  plotOut = plotly::hide_colorbar(p = plotOut)
  plotOut = plotly::hide_legend(p = plotOut)
  plotOut = plotly::config(p = plotOut, displayModeBar=F, editable=T)
  return(plotOut)
}
