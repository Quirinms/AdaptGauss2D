plotDataGridDensity3D = function(Data, XKernel, YKernel,
                               ContinuousDataPDE, EmpiricDataPDE,
                               #Means, Covariances, MainAxesAngle,
                               Colors, Cls, Camera = NULL,
                               ShowScatter = TRUE, #ShowMarkers = NULL,
                               AxNames = c("X", "Y"),
                               Source = "F1", Debug = F){
  # DESCRIPTION
  #
  #
  # INPUT
  # Data[1:n, 1:2]                 Numeric matrix with n observations and 2 features.
  # XKernel[1:x]                   Numeric vector defining domain of x axis.
  # YKernel[1:x]                   Numeric vector defining domain of x axis.
  # ContinuousDataPDE[1:x, 1:x]    Numeric matrix with density estimation of
  #                                Data computed as grid on the domain defined by XKernel
  #                                and YKernel in order to plot a continuous density plot.
  # EmpiricDataPDE[1:n]            Numeric vector with density estimation of Data
  #                                defined for each datapoint within the Data.
  # Colors[1:n]                    Numerical vector of size n containing the
  #                                colors of each observation.
  # Cls[1:n]                       Numerical vector of size n containing the
  #                                classes of each observation.
  # Camera                         List of attributes concerning the camera
  #                                angle for plotly visualizations.
  # ShowScatter                    Boolean indicating if 3D scatter points are
  #                                shown in plot (TRUE) or not (FALSE).
  # AxNames                        Character vector with names for each
  #                                dimension ax of 2D plot.
  # Source                         Character indicating plot source. Important
  #                                attribute for plotly in shiny in order to keep
  #                                control of specific panels.
  # Debug                          Boolean (Default=FALSE). T: Show developer
  #                                information and warnings in terminal.
  #                                F: Show nothing.
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
  if(missing(ContinuousDataPDE)){
    message("Parameter ContinuousDataPDE is missing. Returning.")
    return()
  }else{
    if(!is.matrix(ContinuousDataPDE)){
      message("Parameter ContinuousDataPDE is not of type matrix. Returning.")
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
    if(!is.vector(EmpiricDataPDE)){
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
  #if(is.null(ShowMarkers)){
  #  ShowMarkers = FALSE
  #}
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
  EmpiricDataPDE = EmpiricDataPDE/sum(EmpiricDataPDE)
  colfunc = colorRampPalette(c("white", "blue"))
  #MyColorGradient = c("#FFFFFF", colfunc(10))
  if(requireNamespace("colorRamps"))
	MyColorGradient = colorRamps::matlab.like(10)
  else
	MyColorGradient = c("#FFFFFF", colfunc(10))
	
  MyColorGradient = c(rep("#FFFFFF", 2), MyColorGradient)
  plotOut = plotly::plot_ly(x = XKernel,
                            y = YKernel,
                            z = ContinuousDataPDE,
                            type = "surface",
                            colors = MyColorGradient,
                            alpha = 0.7,
                            source = Source)
  if(ShowScatter){
    plotOut = plotly::add_markers(p = plotOut,
                                  x = Data[,1],
                                  y = Data[,2],
                                  z = EmpiricDataPDE,
                                  marker = list(size = 3,
                                                color = Colors[Cls]))
  }
  plotOut = plotly::event_register(p = plotOut, event = 'plotly_relayout')
  plotOut = plotly::config(p = plotOut, displayModeBar = F)
  plotOut = plotly::hide_colorbar(p = plotOut)
  plotOut = plotly::hide_legend(p = plotOut)
  plotOut = plotly::layout(p = plotOut,
                           title = "3D Empirical PDE",
                           scene = list(xaxis = list(title = AxNames[1], fixedrange =  T),
                                        yaxis = list(title = AxNames[2], fixedrange = T),
                                        zaxis = list(title = "Probability"),
                                        camera = Camera))
  return(plotOut)
}
