plotDataDensity = function(Data, Cls, CurrGauss, Colors,
                           ContinuousDataPDE,
                           Means, Covariances, MainAxesAngle,
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
  # ContinuousDataPDE
  # XKernel[1:x]       Numeric vector defining domain of x axis.
  # YKernel[1:x]       Numeric vector defining domain of x axis.
  # AxNames            Character vector with names for each dimension ax of 2D
  #                    plot.
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
  # DataVisualizations::PmatrixColormap
  #MyColorGradient = colorRamps::matlab.like(5) # do not choose to many colors, so that 0 density values are displayed white
  #MyColorGradient = c("#FFFFFF", "#FFFFFF", MyColorGradient) # White out the first values so that the background appears white
  # Compute the density for a continuous plot of the complete GMM model
  plotOut = plotly::plot_ly(x = XKernel,
                            y = YKernel,
                            type = "contour",
                            colors = AdaptGauss2D::PmatrixColormap,
                            source = Source)
  plotOut = plotly::add_contour(p = plotOut,
                                z = ContinuousDataPDE)
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
  if(ShowScatter){
    plotOut = plotly::add_markers(p = plotOut,
                                  x = Data[,1],
                                  y = Data[,2],
                                  color = Colors[Cls[]],
                                  marker = list(size = 3, color = "black"))#, colors = Colors[1:length(unique(Cls))])
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
                           title  = "2D Data Density",
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
