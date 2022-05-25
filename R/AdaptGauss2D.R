AdaptDunes = function(Data,
                      Means = NULL, CovarianceMatrices = NULL, Weights = NULL,
                      Cls = NULL, Debug = F, dbt = F,
                      WorkingDirectory = getwd()){
  # DESCRIPTION
  # GUI for manually fitting a multivariate (2 dimensional) gaussian mixture
  # to a dataset.
  #
  # USAGE
  # V = AdaptDunes(Data)
  # V = AdaptDunes(Data, Cls = Cls)
  # V = AdaptDunes(Data, Means = Means, CovarianceMatrices = CovarianceMatrices,
  #                Weights = Weights)
  #
  # INPUT
  # Data[1:n, 1:2]        Dataset with two variables for columns that will be
  #                       used for fitting.
  # Means                 List of l numerical 2D vectors.
  # CovarianceMatrices    List of l numerical 2x2 matrices.
  # Weights[1:l]          Numerical vectors of length l.
  # Cls[1:n]              Numerical vectors with class labels for each
  #                       observation.
  # Debug                 Boolean (Default=FALSE). T: Show developer information
  #                       and warnings in terminal. F: Show nothing.
  # dbt                   Boolean: TRUE => databionic team usage and
  #                                        functionality.
  #                                FALSE => external usage with standard
  #                                         utilities.
  #                                Delete all dbt stuff when publishing code.
  # WorkingDirectory      Character indicating working directory for saving
  #                       settings.
  #
  # OUTPUT
  # List of
  #   $Means                     List of l numerical 2D vectors. Each vector
  #                              representing the expectation values for one
  #                              Gaussian component.
  #   $CovarianceMatrices        List of l numerical 2x2 matrices.
  #   $Weights                   Numerical vector of size l containing the
  #                              weight for each Gaussian component within the
  #                              mixture.
  #   $PrincipalComponentAxis    List of l numerical 2D vectors. Each matrix
  #                              representing the prinicpal component axis of
  #                              the covariance matrix of a Gaussian component.
  #   $Angle                     List of l integers. Each integer represents the
  #                              angle of the first principal component ax of
  #                              the covariance matrix.
  #   $Cls                       Numerical vector of size n containing the
  #                              classes of each observation.
  #
  # Author: QS 20.05.2022
  #

  #----------------------------------------------------------------------------#
  #                           Error Capturing
  #----------------------------------------------------------------------------#
  if(missing(Data)){
    message("Parameter Data is missing. Returning.")
    return()
  }else{
    if(!is.matrix(Data)){
      message("Parameter Data is not of type matrix. Returning.")
      return()
    }else if(dim(Data)[2] < 2){
      message("Parameter Data is required to contain two or more features. Returning.")
      return()
    }else if(dim(Data)[1] <= 30){
      message("Parameter Data must contain 30 observations or more. Returning.")
      return()
    }
  }
  if(is.null(colnames(Data))){
    message("No column names in Data were found. Setting the names of the axes
            to ")
    AxNames = c("Dimension X", "Dimension Y")
  }else{
    AxNames = colnames(Data)[1:2]
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
  if(!is.null(CovarianceMatrices)){
    if(!is.list(CovarianceMatrices)){
      message("Parameter CovarianceMatrices is not of type list. Returning.")
      return()
    }else{
      for(i in 1:length(CovarianceMatrices)){
        if(!is.matrix(CovarianceMatrices[[i]])){
          message("Parameter CovarianceMatrices can only contain matrices. Returning.")
          return()
        }else if((dim(CovarianceMatrices[[i]])[1] != 2) | (dim(CovarianceMatrices[[i]])[2] != 2)){
          message("Parameter CovarianceMatrices can only contain matrices of dimension 2x2. Returning.")
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
  if(!is.logical(Debug)){
    message("Parameter Debug is not a logical type. Returning.")
    return()
  }

  gaussColors = c("blue", "red", "darkgoldenrod", "darkgreen", "magenta", "grey",
                  "cyan", "hotpink", "mediumturquoise", "yellow", "dodgerblue", "lightred")
  #gaussColors2 = grDevices::col2rgb(gaussColors, alpha = T)
  #emSelection = c("Baggenstos", "Mixtools")



  #----------------------------------------------------------------------------#
  #                                    UI
  #----------------------------------------------------------------------------#
  ui = fillPage(
    #chooseSliderSkin("Flat", color = "#212729"), #chooseSliderSkin("Modern"),
    #setSliderColor(color = "black", 1:20),
    tags$style("[type = 'number'] {font-size:12px; height:25px;}"),
    # Top panel: "3D Visualizations - Feature 1 - Feature 2 - 2D Visualizations
    wellPanel(style = 'height:90px;',
              fluidRow(
                style = 'padding-left:50px;padding-right:50px;padding-top:5px;padding-bottom:5px;',
                column(4, titlePanel("3D Visualizations")),
                column(2, uiOutput("uiFeature1")),
                column(2, uiOutput("uiFeature2")),
                column(4, titlePanel("2D Visualizations"))
                )),
    # Plot panel: Left plot (3D) - EM Panel / GMM Components - Right plot (2D)
    fluidRow(#style =  'height:40vh',
      style = 'padding-top:5px;padding-bottom:5px;padding-right:5px;',
      column(4, plotlyOutput("leftView",  height = '55vh')
             ),
      column(4,
             fluidRow(
               titlePanel("Select component of GMM"),
               wellPanel(style = 'padding-top:2px;padding-bottom:2px',
                         fluidRow(style="vertical-align:bottom;",#align:top???
                                  column(4,
                                         fluidRow(
                                           style = 'padding-left:20px;padding-right:20px;padding-top:25px;padding-bottom:5px;',
                                           actionButton("remGaussButton", " Remove Current", width = "100%", icon = icon("minus")))
                                         ),
                                  column(4,
                                         fluidRow(
                                           style = 'padding-left:20px;padding-right:20px;padding-top:5px;padding-bottom:5px;',
                                           uiOutput("gaussianNo"),
                                         )
                                         ),
                                  column(4,
                                         fluidRow(
                                           style = 'padding-left:20px;padding-right:20px;padding-top:25px;padding-bottom:5px;',
                                           actionButton("addGaussButton", " Add", width  = "100%", icon = icon("plus")))
                                         )
                                  ),
                         ),
               ),
             uiOutput("control"),
             fluidRow(column(4, uiOutput("controlAngle"), offset = 4))
             ),
      column(4, plotlyOutput("rightView", height  = '55vh')
             )
    ),
    # Manual control panel
    wellPanel(style = 'height:400px;',
              fluidRow(
                #style = 'padding-left:50px;padding-right:50px;padding-top:5px;padding-bottom:5px;',
                column(4,
                       wellPanel(uiOutput("buttonControlLeft")),
                       fluidRow(style = 'padding-left:50px;padding-right:50px;padding-top:5px;padding-bottom:5px;',
                                shiny::column(6,
                                              fluidRow(actionButton("loadCls",
                                                                    "Load Classification",
                                                                    width = "100%", icon = icon("fas fa-file"))),
                                              fluidRow(actionButton("saveAdaptDunes",
                                                                    "Save Current Setting",
                                                                    width = "100%", icon = icon("fas fa-file"))),
                                              fluidRow(actionButton("quitButton",
                                                                    "Quit",
                                                                    width = "100%", icon = icon("door-open")))
                       ))),
                column(4, titlePanel("Expectation Maximization"),
                       wellPanel(style = 'padding-top:5px;padding-bottom:5px',
                                 #conditionalPanel(condition = "dbt==T",
                                 #                  fluidRow(style = 'padding-top:3px;padding-bottom:0px',
                                 #                          column(12, selectInput("EM", "Select EM-Algorithm", choices = emSelection, width = "100%")))),

                                 fluidRow(style = 'padding-top:0px;padding-bottom:2px',
                                          column(4, actionButton("execEM", " EM", width = "100%", icon = icon("calculator"))),
                                          #column(4, actionButton("execEM", " EM", width = "100%", icon = icon("calculator"))),
                                          column(2),
                                          column(3, actionButton("undoEM", "Undo ", width = "100%", icon = icon("undo"))),
                                          column(3, actionButton("redoEM", "Redo ", width = "100%", icon = icon("redo")))
                                 ),
                                 fluidRow(style = 'padding-top:2px;padding-bottom:0px',
                                          column(4, checkboxInput("EMSplitModes", "EM Add/Rem. Modes", value = FALSE, width = "100%")),
                                          column(4, numericInput("EMSteps",
                                                                 "EM Steps",
                                                                 min = 1, max = 200,
                                                                 step = 1,
                                                                 value = 10)),
                                          column(4, textOutput("RMS")),
                                 ),
                                 fluidRow(style = 'padding-top:0px;padding-bottom:0px',
                                          column(6, actionButton("normOth",
                                                                 "Norm Current Weight",
                                                                 width = "100%", icon = icon("balance-scale"))),
                                          column(6, actionButton("normAll",
                                                                 "Norm All Weights",
                                                                 width = "100%", icon = icon("balance-scale"))),
                                 ),
                       )),
                column(4,
                       wellPanel(uiOutput("buttonControlRight")),
                       )
              )
    )
  )

  #----------------------------------------------------------------------------#
  #                                 Server
  #----------------------------------------------------------------------------#
  server<-function(input, output, session){
    if(Debug == F){
      options(warn = -1)    # Do not show any warnings at all if 'Debug' modus is disabled
    }

    #----------------------------------------------------------------------------#
    #                           Initialization
    #----------------------------------------------------------------------------#
    withProgress(message = 'Initializing Adapt Dunes', value = 0, {
      NumStepsInit = 5
      incProgress(1/NumStepsInit, detail = paste("Doing step", 1))
      Features = c()
      if(is.null(colnames(Data))){
        for(i in 1:dim(Data)[2]){
          Features = c(Features, paste0("Feature", i))
        }
      }else{
        Features = colnames(Data)
      }

      Feature1     = Features[1]
      Feature2     = Features[2]
      FeatureIdx1  = 1
      FeatureIdx2  = 2

      if(is.null(Cls)){
        Cls = rep(0, dim(Data)[1])
      }

      # Sample the data (important for Big Data)
      #----------------------------------------------------------------------------#
      # Sample data for EM (use enough to reach an accurate estimate)
      dataSampleIdx = computeSample4EM(nrow(Data), BigDataMaximum = 0.6)
      DataEM        = Data[dataSampleIdx, c(FeatureIdx1, FeatureIdx2)]

      incProgress(1/NumStepsInit, detail = paste("Doing step", 2))

      # Sample data for visualizations (use a minimum of data to reduce
      # computational power but sample smart, so that the structures are visualized
      # proportional)
      dataSampleIdx2 = computeSample4Visualization(nrow(Data))
      DataVis        = Data[dataSampleIdx2, c(FeatureIdx1, FeatureIdx2)]
      incProgress(1/NumStepsInit, detail = paste("Doing step", 3))
      # Density estimation of Data
      #----------------------------------------------------------------------------#
      # Estimate the density both for the data and the domain
      DensityMap        = DataVisualizations::SmoothedDensitiesXY(X = Data[,1], # Uses SDH from [Eilers/Goeman, 2004]
                                                                  Y = Data[,2],
                                                                  nbins = 100)
      ContinuousDataPDE = DensityMap$hist_F_2D/sum(DensityMap$hist_F_2D)        # Density estimation on domain => 3D cont. Density plot
      EmpiricDataPDE    = DensityMap$Densities/sum(DensityMap$Densities)        # Density estimation for data => 3D Scatter plot
      XKernel           = DensityMap$Xkernels
      YKernel           = DensityMap$Ykernels
      #
      #incProgress(1/NumStepsInit, detail = paste("Doing step", 4))
      # PDEScatter?
    })
    #----------------------------------------------------------------------------#
    # Root Mean Squared Deviation
    meanX     = mean(Data[,FeatureIdx1])
    meanY     = mean(Data[,FeatureIdx2])
    DataMean  = c(meanX, meanY)
    SDX       = sd(Data[,FeatureIdx1])
    SDY       = sd(Data[,FeatureIdx2])
    DataSigma = diag(c(SDX, SDY))
    Fi        = mixtools::dmvnorm(y = Data[,c(FeatureIdx1, FeatureIdx2)], mu = DataMean, sigma = DataSigma)
    RMSE0     = sqrt(sum((Fi - EmpiricDataPDE)^2))/100

    # GMM
    #----------------------------------------------------------------------------#
    # If any parameter is not set, estimate a matching number of GMM components
    # and matching settings with EM (here Baggenstoss is used and recommended)
    # Case 1: There are parameters given
    # => Check for correct input and start algorithm
    if(all(!is.null(Means), !(is.null(CovarianceMatrices)), !is.null(Weights))){
      if(!is.list(Means)){
        print("Means is not a list. Returning")
        print("close App")
        stopApp()
      }
      if(!is.list(CovarianceMatrices)){
        print("CovarianceMatrices is not a list. Returning")
        print("close App")
        stopApp()
      }
      if(!is.vector(Weights)){
        print("Weights is not a vector. Returning")
        print("close App")
        stopApp()
      }
      for(i in 1:length(Means)){
        if((length(Means[[i]]) != 2) || (!is.vector(Means[[i]]))){
          print("Means must only contain 2x2 vectors.")
          print("close App")
          stopApp()
        }
      }
      for(i in 1:length(CovarianceMatrices)){
        if((dim(CovarianceMatrices[[i]])[1] != 2) || (dim(CovarianceMatrices[[i]])[1] != dim(CovarianceMatrices[[i]])[2])){
          print("CovarianceMatrices must only contain 2x2 matrices.")
          print("close App")
          stopApp()
        }
      }
      CovMat = CovarianceMatrices
      tmpMainAxesAngles = lapply(CovMat, function(x) covariance2AxesAngle(x))
    }
    # Case 2: Not all parameters are given
    # => Notify user and stop App
    if(any(!is.null(Means), !is.null(CovarianceMatrices), !is.null(Weights))){
      print("Set all parameters or none.")
      print("close App")
      stopApp()
    }
    # Case 3: No parameters are given
    # => Continue with results from Expectation Maximization from Baggenstoss
    if(any(is.null(Means), is.null(CovarianceMatrices), is.null(Weights))){
      V = EstimateNumberOfModes(Data = Data[,1:2], MaxModeNo = 7)
      tmpRes = ad_interfaceEM(Data = DataEM, Mu = NULL, Sigma = NULL, Lambda = NULL,
                              Nmodes = V$ModesNo, Addmodes = 0, Nit = 50, Verbose = 0)
      Means             = tmpRes$Mu
      CovMat            = tmpRes$Sigma
      Weights           = tmpRes$Lambda
      tmpMainAxesAngles = tmpRes$MainAxesAngle
    }

    #----------------------------------------------------------------------------#
    if(Debug){
      cat(file=stderr(), paste0("Total no. cases: ",         dim(Data)[1]) ,"\n")
      cat(file=stderr(), paste0("Total no. features: ",      dim(Data)[2]) ,"\n")
      cat(file=stderr(), paste0("Name of feature 1: ",       Feature1) ,"\n")
      cat(file=stderr(), paste0("Name of feature 2: ",       Feature2) ,"\n")
      cat(file=stderr(), paste0("Sum of ContinuousDataPDE ", sum(ContinuousDataPDE)) ,"\n")
      cat(file=stderr(), paste0("Sum of EmpiricDataPDE ",    sum(EmpiricDataPDE)) ,"\n")
      cat(file=stderr(), paste0("Number of gaussians: ",     length(Means)) ,"\n")
      cat(file=stderr(), paste0("Sum of weights: ",          sum(Weights)) ,"\n")
    }

    #--------------------------------------------------------------------------#
    #          Reactive Values
    #--------------------------------------------------------------------------#
    # Reactive values for interactive functionality
    values                      <<- reactiveValues()
    #values$errorMessage         = ""
    values$updateEllipsoid      = 0
    values$lastKnownShapeStart  = matrix(0,ncol=2, nrow = 15)
    values$shapes               = list()
    values$shapeText            = c()
    values$Annotation           = rbind()
    values$dataScatter          = TRUE
    values$DotOrGrid            = "Grid Density"
    values$leftPlot             = "Data Estimation"
    values$leftPlot2            = "Orig. Classes"
    values$lastPlotExecution    = Sys.time() - 501
    values$ShowMarkersLeft      = NULL
    values$rightPlot            = 1
    values$showAxis             = "No axis"
    values$showEllipsoids       = "No ellipsoids"
    values$showScatter          = "No scatter"
    values$rightCls             = "Orig. Classes"

    # Reactive values for internal updates
    befehl                        <<- reactiveValues()
    befehl$updateGaussianNo       =   0
    befehl$updateCurrGauss        =   0
    befehl$drawUIControl          =   0
    befehl$drawGaussianComponents =   0
    #befehl$updateSlider           =   0
    befehl$updateNumericInputs    =   0
    befehl$plot                   =   0
    befehl$updateVisualizations   =   0
    befehl$updateClassification   =   0
    befehl$updateEllipsoid        =   0
    befehl$renderRMS              =   0
    befehl$buttonClasses          =   0

    #--------------------------------------------------------------------------#
    #          Global variables
    #--------------------------------------------------------------------------#
    #
    # User defined inputs
    OriginalCls          = Cls
    #
    # Variables based on Data
    MinDataVal           = round(min(min(Data[,1]), min(Data[,2])), 2)    # Minimum for diverse (e.g., plots are all scaled the same)
    MaxDataVal           = round(max(max(Data[,1]), max(Data[,2])), 2)    # Maximum for diverse (e.g., range of parameters Mean)
    LimitMean            = c(MinDataVal, MaxDataVal)                      # Limit of Mean
    LimitAngle           = c(0, 359.99)                                   # Limit of Angle for the 1st PCA Axis
    LimitMainAxes        = c(0, 0.7*MaxDataVal)                           # Limit for size of PCA axes
    #
    # Variables based on GMM for shiny application
    CurrGauss            = 1                                              # Current Gauss to work with (mainly for showing in app)
    NumGauss             = length(Means)                                  # Total number of Gaussian components
    #
    # Variables based on GMM
    if(NumGauss > length(gaussColors)){ # If there are more components than colors, cut them (current length is 12!)
      NumGauss      = length(gaussColors)
      Means         = Means[1:NumGauss]                                   # List of Means
      CovMats       = CovMat[1:NumGauss]                                  # List of Covariance Matrices
      Weights       = Weights[1:NumGauss]                                 # Vector of the Weights
      MainAxesAngle = tmpMainAxesAngles                                   # List of Princ. Comp. Axes + Angle (of 1st PCA Axis)
    }else{
      CovMats       = CovMat                                              # List of Covariance Matrices
      MainAxesAngle = tmpMainAxesAngles                                   # List of Princ. Comp. Axes + Angle (of 1st PCA Axis)
    }
    #
    # Root Mean Square Deviation (Normalization with RMSE of one gaussian model)
    RMS                  = 0
    #
    # Variables based on GMM and Data
    V                    = computeGMMClassification(DataEM, Means,
                                                    CovMats, Weights)
    Classification       = V$Classification                               # Classification of Data by current GMM
    CurrDataDensity      = V$CurrDataDensity                              # rep(0, nrow(Data))# density for each datapoint within the mixture
    #
    # Variables for EM
    EMNumIterations      = 10
    EMSplitModes         = 0
    #
    # Variables for Visualizations
    GridDensity          = computeGridDensity(XKernel, YKernel,
                                              Means, CovMats, Weights)
    DotDensity           = computeDotDensity(DataVis, Means, CovMats, Weights)
    #
    # Shiny settings for plots
    cameraLeft           = NULL
    cameraRight          = NULL
    # Backup variables
    BackUpMeans          = list()                                         # Backup for undoing results of the EM
    BackUpCovMats        = list()                                         #
    BackUpMainAxesAngle  = list()                                         #
    BackUpWeights        = list()                                         #
    BackUpCurrGauss      = c()                                            #
    CurrBackUpPosition   = 0
    # Shiny control variables
    iBefehl              = 0                                              #
    MonitorStopReactions = F                                              # Boolean value to deploy a reactive Timer
    RightPlotControlStopReactions = F
    ControlPlotReaction1 = F
    ControlPlotReaction2 = F
    ControlPlotReaction3 = F
    ControlEM            = F
    ControlAddGauss      = F
    ControlRemoveGauss   = F
    ControlUndoRedoEM    = F
    ControlNormAll       = F
    ControlNormOthers    = F
    Control6             = F
    ControlSI            = F

    #--------------------------------------------------------------------------#
    #     Dynamically rendered components of UI
    #--------------------------------------------------------------------------#
    # => Changing Action Buttons
    # => Gaussian Components Selection as Slider
    # => Slider Inputs (Weights, Means, Main Axes, Angle of first main ax)
    # => Numeric Inputs (Weights, Means, Main Axes, Angle of first main ax)

    output$uiFeature1 = renderUI({
      uiFeature1Selection = selectInput(inputId   = "Feature1",
                                        label     = "Feature 1",
                                        choices   = Features,
                                        selected  = Feature1,
                                        multiple  = FALSE,
                                        selectize = TRUE,
                                        width     = "100%")
      uiFeature1Selection
    })

    output$uiFeature2 = renderUI({
      uiFeature2Selection = selectInput(inputId   = "Feature2",
                                        label     = "Feature 2",
                                        choices   = Features,
                                        selected  = Feature2,
                                        multiple  = FALSE,
                                        selectize = TRUE,
                                        width     = "100%")
      uiFeature2Selection
    })

    output$gaussianNo = renderUI({
      befehl$drawGaussianComponents
      GaussianNoRow =
        fluidRow(style="vertical-align:bottom;",
                 shiny::numericInput("gaussianNo",
                                     paste0("Gaussian No.: ", CurrGauss, " of ", NumGauss),
                                     min = 1, max = NumGauss,
                                     value = CurrGauss,
                                     step = 1,
                                     width="100%")
                 #shiny::sliderInput("gaussianNo",           # Slider Input Gaussian (component) Number
                #                    paste0("Gaussian No.: ", CurrGauss),
                #                    min = 1, max = NumGauss,
                #                    value = CurrGauss,
                #                    step = 1,
                #                    width="100%")
        )
      GaussianNoRow
    })

    output$buttonControlLeft = renderUI({
      buttonControlLeft =
        fluidRow(style="vertical-align:bottom;", #align:top???
                 column(6, # conditionalPanel
                        fluidRow(actionButton("changeLeftPlot", paste0(values$leftPlot, " (P)"), value = T, width = "100%")),
                        fluidRow(uiOutput("classButtonLeft")),
                        #fluidRow(checkboxInput("showMarkersLeft", "Show Markers", value = F, width = "100%"))
                 ),
                 column(6,
                        fluidRow(actionButton("ChangeDotOrGrid", paste0(values$DotOrGrid, " (P)"), value = T, width = "100%")),
                        fluidRow(actionButton("dataScatter", "Data Scatter", value = TRUE, width = "100%"))
                 ),
        )
      buttonControlLeft
    })

    output$classButtonLeft <- renderUI({
      befehl$buttonClasses
      if(!all(OriginalCls == rep(0, dim(Data)[1]))){
        actionButton("changeLeftPlot2", values$leftPlot2, value = T, width = "100%")
      }
    })

    output$buttonControlRight = renderUI({
      buttonControlRight =
        fluidRow(style="vertical-align:bottom;",#align:top???
                 column(3,
                        fluidRow(actionButton("changeRightPlot1", "Scatter (P)",        value = T, width = "100%")),
                        fluidRow(actionButton("showAxis",       values$showAxis,       value = T, width = "100%"))
                 ),
                 column(3,
                        fluidRow(actionButton("changeRightPlot2", "Model Density (P)",  value = T, width = "100%")),
                        fluidRow(actionButton("showEllipsoids", values$showEllipsoids, value = T, width = "100%"))
                 ),
                 column(3,
                        fluidRow(actionButton("changeRightPlot3", "Data Density (P)",   value = T, width = "100%")),
                        fluidRow(actionButton("showScatter",    values$showScatter,    value = T, width = "100%"))
                 ),
                 column(3,
                        fluidRow(actionButton("changeRightPlot4", "Classification (P)", value = T, width = "100%")),
                        fluidRow(uiOutput("classButtonRight"))
                 )
        )
      buttonControlRight
    })

    output$classButtonRight <- renderUI({
      befehl$buttonClasses
      if(!all(OriginalCls == rep(0, dim(Data)[1]))){
        actionButton("changeRightCls", values$rightCls, value = T, width = "100%")
      }
    })


    output$control = renderUI({
      befehl$drawUIControl
      #gIds = paste0("gaussWeight", 1:NumGauss)
      LowerControl =
        shiny::fluidRow(
          titlePanel(paste0("Current component: ", CurrGauss, " (", gaussColors[CurrGauss], ")")),
          wellPanel(
            shiny::fluidRow(
              shiny::column(4,
                            shiny::numericInput("NumericMean1",  # Numeric Input Center X
                                                paste0("Mean of ", AxNames[1]),
                                                min = LimitMean[1], max = LimitMean[2],
                                                step = 0.01,
                                                value = Means[[CurrGauss]][1]),
                            style=paste0("color:",gaussColors[CurrGauss])
                            ),
              shiny::column(4,
                            shiny::numericInput("NumericMainAxis1", # Numeric Input PCA 1
                                                "Main Axis 1",
                                                min = LimitMainAxes[1], max = LimitMainAxes[2],
                                                step = 0.01,
                                                value = MainAxesAngle[[CurrGauss]][1]),
                            style=paste0("color:",gaussColors[CurrGauss])
                            ),
              shiny::column(4,
                            shiny::numericInput("NumericWeight",        # Numeric Input Weight
                                                paste0("Weight: ", CurrGauss),
                                                min = 0, max = 1,
                                                step = 0.01,
                                                value = Weights[CurrGauss]),
                            style=paste0("color:",gaussColors[CurrGauss])
                            ),
            ),
          shiny::fluidRow(
            shiny::column(4,
                          shiny::numericInput("NumericMean2",    # Numeric Input Center Y
                                              paste0("Mean of ", AxNames[2]),
                                              min = LimitMean[1], max = LimitMean[2],
                                              step = 0.01,
                                              value = Means[[CurrGauss]][2]),
                          style=paste0("color:",gaussColors[CurrGauss])
                          ),
            shiny::column(4,
                          shiny::numericInput("NumericMainAxis2",    # Numeric Input PCA 2
                                              "Main Axis 2",
                                              min = LimitMainAxes[1], max = LimitMainAxes[2],
                                              step = 0.01,
                                              value = MainAxesAngle[[CurrGauss]][2]),
                          style=paste0("color:",gaussColors[CurrGauss])
                          ),
            shiny::column(4,
                          shiny::numericInput("NumericAngle",   # Numeric Input Angle
                                              "Angle",
                                              min = LimitAngle[1], max = LimitAngle[2],
                                              step = 0.01,
                                              value = decimal2CompassDegree(MainAxesAngle[[CurrGauss]][3])),
                          style=paste0("color:",gaussColors[CurrGauss])
                          )
            )
          )
        )
      LowerControl
    })

    output$controlAngle = renderUI({
      befehl$drawUIControl
      #gIds = paste0("gaussWeight", 1:NumGauss)
      LowerControl = shiny::column(12,
                                   shinyWidgets::knobInput(
                                     inputId         = "SliderAngle",
                                      label           = "Angle of Main Axis 1",
                                      value           = decimal2CompassDegree(MainAxesAngle[[CurrGauss]][3]),
                                      min             = LimitAngle[1], max = LimitAngle[2],
                                      displayPrevious = TRUE,
                                      lineCap         = "round",
                                      step            = 0.1,
                                      width           = "100%",
                                      fgColor         = "steelblue",
                                      inputColor      = "steelblue"
                                    ),
                                    style=paste0("color:",gaussColors[CurrGauss]))
      LowerControl
    })

    #--------------------------------------------------------------------------#
    #                Internal functions
    #--------------------------------------------------------------------------#

    create_backup = function(){
      #print(paste0("create_backup", tmpOld))
      if(length(BackUpMeans[[CurrBackUpPosition]]) == length(Means)){    # Check if current setting unequals new setting, implying need for backup
        NoBackUpRequired = TRUE
        for(i in 1:length(BackUpMeans[[CurrBackUpPosition]])){
          A = all(round(BackUpMeans[[CurrBackUpPosition]][[i]],2) == round(Means[[i]], 2))
          B = all(round(BackUpCovMats[[CurrBackUpPosition]][[i]],1) == round(CovMats[[i]],1))
          C = all(round(BackUpMainAxesAngle[[CurrBackUpPosition]][[i]][1:2],2) == round(MainAxesAngle[[i]][1:2],2))
          C2 = all(round(BackUpMainAxesAngle[[CurrBackUpPosition]][[i]][3],0) == round(MainAxesAngle[[i]][3],0))
          D = (round(BackUpWeights[[CurrBackUpPosition]][i],2) == round(Weights[i],2))
          NoBackUpRequired = NoBackUpRequired & A & B & C & D
        }
        if(NoBackUpRequired == FALSE){
          purge_future_history()         # Remove all Backups which were ahead of the current setting, since another path was chosen
          execute_backup()               # Create backup after the current position
        }
      }else{
        purge_future_history()         # Remove all Backups which were ahead of the current setting, since another path was chosen
        execute_backup()               # Create backup after the current position
      }
    }

    purge_future_history = function(){
      if(length(BackUpMeans) > CurrBackUpPosition){
        intL = length(BackUpMeans) - CurrBackUpPosition
        BackUpMeans         <<- BackUpMeans[1:CurrBackUpPosition]
        BackUpCovMats       <<- BackUpCovMats[1:CurrBackUpPosition]
        BackUpMainAxesAngle <<- BackUpMainAxesAngle[1:CurrBackUpPosition]
        BackUpWeights       <<- BackUpWeights[1:CurrBackUpPosition]
        BackUpCurrGauss     <<- BackUpCurrGauss[1:CurrBackUpPosition]
      }
    }

    execute_backup = function(){
      CurrBackUpPosition                        <<- CurrBackUpPosition + 1
      #tmpNew                        =   length(BackUpMeans) + 1
      #print(paste0("execute_backup", tmpNew))
      BackUpMeans[[CurrBackUpPosition]]         <<- Means
      BackUpCovMats[[CurrBackUpPosition]]       <<- CovMats
      BackUpMainAxesAngle[[CurrBackUpPosition]] <<- MainAxesAngle
      BackUpWeights[[CurrBackUpPosition]]       <<- Weights
      BackUpCurrGauss                           <<- c(BackUpCurrGauss, CurrGauss)
    }

    redo_next_backup = function(){
      if(length(BackUpMeans) >= (CurrBackUpPosition + 1)){
        CurrBackUpPosition <<- CurrBackUpPosition + 1
        CurrGauss          <<- BackUpCurrGauss[CurrBackUpPosition]
        Means              <<- BackUpMeans[[CurrBackUpPosition]]
        CovMats            <<- BackUpCovMats[[CurrBackUpPosition]]
        MainAxesAngle      <<- BackUpMainAxesAngle[[CurrBackUpPosition]]
        Weights            <<- BackUpWeights[[CurrBackUpPosition]]
        NumGauss           <<- length(Means)
        iBefehl                       <<- iBefehl+1
        befehl$updateCurrGauss        <-  iBefehl
        befehl$drawGaussianComponents <-  iBefehl
        befehl$updateClassification   <-  iBefehl
        befehl$updateVisualizations   <-  iBefehl
        #befehl$updateSlider           <-  iBefehl
        befehl$updateNumericInputs    <-  iBefehl
        befehl$updateEllipsoid        <-  iBefehl
        befehl$plot                   <-  iBefehl
        befehl$drawUIControl          <-  iBefehl
      }else{                                                              # Else, there is no backup, notify user and continue
        cat(file = stderr(), "Redo not possible. History is at the first position.\n")
      }
    }

    recover_latest_backup = function(){
      #GoTo       = length(BackUpMeans) - 1
      #NumBackUps = length(BackUpMeans)
      GoTo       = CurrBackUpPosition - 1
      NumBackUps = CurrBackUpPosition
      #print(paste0("recover_latest_backup", NumBackUps))
      if((NumBackUps > 1) & (GoTo > 0)){                                  # If there is a backup, use it
        CurrBackUpPosition <<- CurrBackUpPosition - 1
        CurrGauss     <<- BackUpCurrGauss[CurrBackUpPosition]
        Means         <<- BackUpMeans[[CurrBackUpPosition]]                             # Use Backup to load previous results
        CovMats       <<- BackUpCovMats[[CurrBackUpPosition]]
        MainAxesAngle <<- BackUpMainAxesAngle[[CurrBackUpPosition]]
        Weights       <<- BackUpWeights[[CurrBackUpPosition]]
        NumGauss      <<- length(Means)
        #if(NumBackUps > 1){                                               # Never delete the first entry (always being able to restart)
        #  BackUpCurrGauss     <<- BackUpCurrGauss[1:GoTo]
        #  BackUpMeans         <<- BackUpMeans[1:GoTo]                     # Remove latest result from Backup
        #  BackUpCovMats       <<- BackUpCovMats[1:GoTo]
        #  BackUpMainAxesAngle <<- BackUpMainAxesAngle[1:GoTo]
        #  BackUpWeights       <<- BackUpWeights[1:GoTo]
        #}
        iBefehl                       <<- iBefehl+1
        befehl$updateCurrGauss        <-  iBefehl
        befehl$drawGaussianComponents <-  iBefehl
        befehl$updateClassification   <-  iBefehl
        befehl$updateVisualizations   <-  iBefehl
        #befehl$updateSlider           <-  iBefehl
        befehl$updateNumericInputs    <-  iBefehl
        befehl$updateEllipsoid        <-  iBefehl
        befehl$plot                   <-  iBefehl
        befehl$drawUIControl          <-  iBefehl
      }else{                                                              # Else, there is no backup, notify user and continue
        cat(file = stderr(), "Undo not possible. History is at the first position.\n")
      }
    }

    #--------------------------------------------------------------------------#
    # Last steps before the UI starts and after the init. of all vars and funs
    #--------------------------------------------------------------------------#
    execute_backup()                # Save the start setting

    #--------------------------------------------------------------------------#
    #     Observe update for slider and numeric input
    #--------------------------------------------------------------------------#
    # Every 500ms an update is allowed
    # Delaying time at this point avoids infinite update loops
    validationTimer1 <- reactiveTimer(500)
    validationTimer2 <- reactiveTimer(1000)
    validationTimer3 <- reactiveTimer(1000)
    validationTimer4 <- reactiveTimer(1000)
    validationTimerEM <- reactiveTimer(5000)
    validationTimer5 <- reactiveTimer(5000)
    validationTimer6 <- reactiveTimer(5000)
    validationTimer7 <- reactiveTimer(5000)
    validationTimer8 <- reactiveTimer(5000)
    validationTimer9 <- reactiveTimer(5000)
    validationTimer10 <- reactiveTimer(5000)
    validationTimerC1 <- reactiveTimer(5000)
    validationTimerC2 <- reactiveTimer(5000)
    validationTimerC3 <- reactiveTimer(5000)
    validationTimerC4 <- reactiveTimer(5000)
    validationTimerC5 <- reactiveTimer(5000)
    validationTimerC6 <- reactiveTimer(5000)
    validationTimerSI <- reactiveTimer(5000)
    validationTimerRightPlotControl <- reactiveTimer(2000)



    observe({
      validationTimerRightPlotControl()
      RightPlotControlStopReactions <<- F
    })

    observe({
      validationTimer1()
      MonitorStopReactions <<- F
    })

    observe({
      validationTimer2()
      ControlPlotReaction1 <<- F
    })

    observe({
      validationTimer3()
      ControlPlotReaction2 <<- F
    })

    observe({
      validationTimer4()
      ControlPlotReaction3 <<- F
    })

    observe({
      validationTimerEM()
      ControlEM <<- F
    })

    observe({
      validationTimer5()
      ControlAddGauss <<- F
    })

    observe({
      validationTimer6()
      ControlRemoveGauss <<- F
    })

    observe({
      validationTimer7()
      ControlUndoRedoEM <<- F
    })

    observe({
      validationTimer8()
      ControlNormAll <<- F
    })

    observe({
      validationTimer9()
      ControlNormOthers <<- F
    })

    observe({
      validationTimer10()
      Control6 <<- F
    })

    observe({
      validationTimerSI()
      ControlSI <<- F
    })

    #--------------------------------------------------------------------------#
    # Change features (consider each dimension alone standing)
    #--------------------------------------------------------------------------#
    # If a feature is exchanged for a new one, a complete re-computation of all
    # initiating steps (as if the app would start from scratch) is necessary.
    observe({
      input$Feature1
      if(!is.null(input$Feature1)){
        if(input$Feature1 != Feature1){
          Feature1    <<- as.character(input$Feature1)
          FeatureIdx1 <<- which(Features == Feature1)[1]
          AxNames[1]  <<- Feature1
          DataEM      <<- Data[dataSampleIdx, c(FeatureIdx1, FeatureIdx2)]
          DataVis     <<- Data[dataSampleIdx2, c(FeatureIdx1, FeatureIdx2)]

          # Root Mean Squared Deviation
          meanX     = mean(Data[,FeatureIdx1])
          meanY     = mean(Data[,FeatureIdx2])
          DataMean  = c(meanX, meanY)
          SDX       = sd(Data[,FeatureIdx1])
          SDY       = sd(Data[,FeatureIdx2])
          DataSigma = diag(c(SDX, SDY))
          Fi        = mixtools::dmvnorm(y = Data[,c(FeatureIdx1, FeatureIdx2)], mu = DataMean, sigma = DataSigma)
          RMSE0     = sqrt(sum((Fi - EmpiricDataPDE)^2))

          BackUpMeans          <<- list()                                         # Backup for undoing results of the EM
          BackUpCovMats        <<- list()                                         #
          BackUpMainAxesAngle  <<- list()                                         #
          BackUpWeights        <<- list()                                         #
          BackUpCurrGauss      <<- c()                                            #
          CurrBackUpPosition   <<- 0

          # Variables based on Data
          MinDataVal           <<- round(min(min(Data[,FeatureIdx1]), min(Data[,FeatureIdx2])), 2)    # Minimum for diverse (e.g., plots are all scaled the same)
          MaxDataVal           <<- round(max(max(Data[,FeatureIdx1]), max(Data[,FeatureIdx2])), 2)    # Maximum for diverse (e.g., range of parameters Mean)
          LimitMean            <<- c(MinDataVal, MaxDataVal)

          # Call EM
          V   = EstimateNumberOfModes(Data = Data[,c(FeatureIdx1, FeatureIdx2)], MaxModeNo = 7)
          res = ad_interfaceEM(Data = DataEM, Mu = NULL, Sigma = NULL,
                               Lambda = NULL, Nmodes = V$ModesNo, Addmodes = 0, Nit = 50)

          # Variables based on GMM
          Means         <<- res$Mu
          CovMats       <<- res$Sigma
          Weights       <<- res$Lambda
          MainAxesAngle <<- res$MainAxesAngle
          NumGauss      <<- length(Means)
          CurrGauss     <<- 1

          execute_backup()

          DensityMap        <<- DataVisualizations::SmoothedDensitiesXY(X = Data[,FeatureIdx1], # Uses SDH from [Eilers/Goeman, 2004]
                                                                        Y = Data[,FeatureIdx2],
                                                                        nbins = 100)
          ContinuousDataPDE <<- DensityMap$hist_F_2D/sum(DensityMap$hist_F_2D)        # Density estimation on domain => 3D cont. Density plot
          EmpiricDataPDE    <<- DensityMap$Densities/sum(DensityMap$Densities)        # Density estimation for data => 3D Scatter plot
          XKernel           <<- DensityMap$Xkernels
          YKernel           <<- DensityMap$Ykernels

          iBefehl                       <<- iBefehl+1
          befehl$updateClassification   <-  iBefehl
          befehl$updateVisualizations   <-  iBefehl
          #befehl$updateSlider           <-  iBefehl
          befehl$updateNumericInputs    <-  iBefehl
          befehl$updateCurrGauss        <-  iBefehl
          befehl$updateEllipsoid        <-  iBefehl
          befehl$drawUIControl          <-  iBefehl
          befehl$drawGaussianComponents <-  iBefehl
        }
      }
    })

    observe({
      input$Feature2
      if(!is.null(input$Feature2)){
        if(input$Feature2 != Feature2){
          Feature2    <<- as.character(input$Feature2)
          FeatureIdx2 <<- which(Features == Feature2)[1]
          AxNames[2]  <<- Feature2
          DataEM[,2]  <<- Data[dataSampleIdx, FeatureIdx2]
          DataVis[,2] <<- Data[dataSampleIdx2, FeatureIdx2]

          # Root Mean Squared Deviation
          meanX     = mean(Data[,FeatureIdx1])
          meanY     = mean(Data[,FeatureIdx2])
          DataMean  = c(meanX, meanY)
          SDX       = sd(Data[,FeatureIdx1])
          SDY       = sd(Data[,FeatureIdx2])
          DataSigma = diag(c(SDX, SDY))
          Fi        = mixtools::dmvnorm(y = Data[,c(FeatureIdx1, FeatureIdx2)], mu = DataMean, sigma = DataSigma)
          RMSE0     = sqrt(sum((Fi - EmpiricDataPDE)^2))

          BackUpMeans          <<- list()                                         # Backup for undoing results of the EM
          BackUpCovMats        <<- list()                                         #
          BackUpMainAxesAngle  <<- list()                                         #
          BackUpWeights        <<- list()                                         #
          BackUpCurrGauss      <<- c()                                            #
          CurrBackUpPosition   <<- 0

          # Variables based on Data
          MinDataVal           <<- round(min(min(Data[,FeatureIdx1]), min(Data[,FeatureIdx2])), 2)    # Minimum for diverse (e.g., plots are all scaled the same)
          MaxDataVal           <<- round(max(max(Data[,FeatureIdx1]), max(Data[,FeatureIdx2])), 2)    # Maximum for diverse (e.g., range of parameters Mean)
          LimitMean            <<- c(MinDataVal, MaxDataVal)

          # Call EM
          V   = EstimateNumberOfModes(Data = Data[,c(FeatureIdx1, FeatureIdx2)], MaxModeNo = 7)
          res = ad_interfaceEM(Data = DataEM, Mu = NULL, Sigma = NULL,
                               Lambda = NULL, Nmodes = V$ModesNo, Addmodes = 0, Nit = 50)

          # Variables based on GMM
          Means         <<- res$Mu
          CovMats       <<- res$Sigma
          Weights       <<- res$Lambda
          MainAxesAngle <<- res$MainAxesAngle
          NumGauss      <<- length(Means)
          CurrGauss     <<- 1

          execute_backup()

          DensityMap        <<- DataVisualizations::SmoothedDensitiesXY(X = Data[,FeatureIdx1], # Uses SDH from [Eilers/Goeman, 2004]
                                                                        Y = Data[,FeatureIdx2],
                                                                        nbins = 100)
          ContinuousDataPDE <<- DensityMap$hist_F_2D/sum(DensityMap$hist_F_2D)        # Density estimation on domain => 3D cont. Density plot
          EmpiricDataPDE    <<- DensityMap$Densities/sum(DensityMap$Densities)        # Density estimation for data => 3D Scatter plot
          XKernel           <<- DensityMap$Xkernels
          YKernel           <<- DensityMap$Ykernels

          iBefehl                       <<- iBefehl+1
          befehl$updateClassification   <-  iBefehl
          befehl$updateVisualizations   <-  iBefehl
          #befehl$updateSlider           <-  iBefehl
          befehl$updateNumericInputs    <-  iBefehl
          befehl$updateCurrGauss        <-  iBefehl
          befehl$updateEllipsoid        <-  iBefehl
          befehl$drawUIControl          <-  iBefehl
          befehl$drawGaussianComponents <-  iBefehl
        }
      }
    })

    # Update the visualizations
    observe({
      befehl$updateVisualizations
      GridDensity <<- computeGridDensity(XKernel, YKernel,
                                         Means, CovMats, Weights)
      DotDensity <<- computeDotDensity(DataVis, Means, CovMats, Weights)
    })

    # Update Numeric Inputs (e.g., gaussian component was changed)
    observe({
      #print("AdaptDunes: Update Slider for M, S and W")
      befehl$updateNumericInputs       #
      befehl$updateCurrGauss    # Gaussian component was updated
      disable("execEM")
      #print("Update Numeric Input M/S/W")
      updateNumericInput(session, "NumericWeight",    value = round(Weights[CurrGauss], 4))
      updateNumericInput(session, "NumericMean1",     value = round(Means[[CurrGauss]][1], 4))
      updateNumericInput(session, "NumericMean2",     value = round(Means[[CurrGauss]][2], 4))
      updateNumericInput(session, "NumericAngle",     value = round(decimal2CompassDegree(MainAxesAngle[[CurrGauss]][3]), 4))
      updateNumericInput(session, "NumericMainAxis1", value = round(MainAxesAngle[[CurrGauss]][1], 4))
      updateNumericInput(session, "NumericMainAxis2", value = round(MainAxesAngle[[CurrGauss]][2], 4))
    })

    # Update the number of iterations in each EM computation
    observe({
      input$EMSteps
      if(is.numeric(input$EMSteps))
        EMNumIterations <<- input$EMSteps
      #updateNumericInput(session, "EMNumIterations", value = EMNumIterations)
    })

    # Update checkbox for EM setting: Allow EM to split modes (TRUE/Check box) or not (FALSE/Empty box)
    observe({
      input$EMSplitModes
      EMSplitModes <<- as.numeric(input$EMSplitModes)
    })

    # Update Numeric Inputs if Slider Inputs change
    observe({
      if(!is.null(input$SliderAngle)){
        updateNumericInput(session, "NumericAngle", value=as.numeric(input$SliderAngle))
      }
    })


    # Update Slider Inputs if Numeric Inputs change
    observe({
      tmpVar = input$NumericAngle
      # Do not update numeric input too fast otherwise infinite loop with mutual updates with slider
      if(MonitorStopReactions==F){
        MonitorStopReactions<<-T
        ControlSI<<-T
        if(!is.null(tmpVar)){
          if(is.numeric(tmpVar)){
            if((as.numeric(tmpVar) >= LimitAngle[1]) & (as.numeric(tmpVar) < LimitAngle[2])){
              updateSliderInput(session, "SliderAngle", value=tmpVar)
            }
          }
        }
      }
    })

    # Update Slider of gaussian components
    observe({
      befehl$updateGaussianNo
      Control6 <<- T
      updateSliderInput(session, "gaussianNo", value = CurrGauss)
    })

    # calculate the current classification
    observe({
      befehl$updateClassification
      disable("execEM")
      # Nota:
      # Prevent user input if computation is taking place
      # (currently necessary (apparently) only here)
      # Extend to other complex computations if necessary
      if(Debug){
        cat(file = stderr(), "Calculating the density and classification\n")
      }
      # Compute classification only on the sample which needs to be visualized (full classification later)
      tmpV      = computeGMMClassification(DataVis, Means, CovMats, Weights)
      tmpCls                 = Classification
      tmpCls[dataSampleIdx2] = tmpV$Classification
      Classification  <<- tmpCls
      #Classification <<- tmpV$Classification
      CurrDataDensity <<- tmpV$CurrDataDensity
      RMS             <<- computeRMS(Data[,c(FeatureIdx1, FeatureIdx2)], Means, CovMats, Weights, EmpiricDataPDE)
      iBefehl          <<- iBefehl+1
      befehl$renderRMS <-  iBefehl
    })

    observe({
      befehl$renderRMS
      output$RMS = renderText({paste0("RMSD%: ", as.character(round(RMS/RMSE0, 4)))})
    })

    # Refreshes values for Weights, Means, Axes and Angles if the Slider is changed
    observe({
      Control6<<-T
      tmpVar = input$NumericWeight
      #print("Refresh Values for Weights")
      #print("AdaptDunes:Refresh Weights")
      if(!is.null(tmpVar)){
        if(is.numeric(tmpVar)){
          if((tmpVar > 0) & (tmpVar < 1)){
            Weights[CurrGauss]            <<- tmpVar
            create_backup()
            iBefehl                       <<- iBefehl+1
            befehl$updateClassification   <-  iBefehl
            befehl$updateVisualizations   <-  iBefehl
            befehl$updateEllipsoid        <-  iBefehl
            befehl$plot                   <-  iBefehl
          }
        }
      }
    })

    observe({
      Control6<<-T
      tmpVar = input$NumericMean1
      #print("Refresh Values for Means")
      #print("AdaptDunes:Refresh Means")
      if(!is.null(tmpVar)){
        if(is.numeric(tmpVar)){
          if((tmpVar > LimitMean[1]) & (tmpVar < LimitMean[2])){
            Means[[CurrGauss]][1]         <<- tmpVar
            create_backup()
            iBefehl                       <<- iBefehl+1
            befehl$updateClassification   <-  iBefehl
            befehl$updateVisualizations   <-  iBefehl
            befehl$updateEllipsoid        <-  iBefehl
            befehl$plot                   <-  iBefehl
          }
        }
      }
    })

    observe({
      Control6 <<- T
      tmpVar = input$NumericMean2
      #print("Refresh Values for Means")
      #print("AdaptDunes:Refresh Means")
      if(!is.null(tmpVar)){
        if(is.numeric(tmpVar)){
          if((tmpVar > LimitMean[1]) & (tmpVar < LimitMean[2])){
            Means[[CurrGauss]][2]         <<- tmpVar
            create_backup()
            iBefehl                       <<- iBefehl+1
            befehl$updateClassification   <-  iBefehl
            befehl$updateVisualizations   <-  iBefehl
            befehl$updateEllipsoid        <-  iBefehl
            befehl$plot                   <-  iBefehl
          }
        }
      }
    })

    observe({
      tmpVar = input$NumericAngle
      #print("Refresh Values for Angle")
      #print("AdaptDunes:Refresh Angle")
      if(!is.null(tmpVar)){
        if(is.numeric(tmpVar)){
          if((tmpVar >= LimitAngle[1]) & (tmpVar <= LimitAngle[2])){
            MainAxesAngle[[CurrGauss]][3] <<- compass2DecimalDegree(as.numeric(tmpVar))
            CovMats <<- lapply(MainAxesAngle,                                      # Update the Covariance Matrix
                               function(x) axesAngle2Covariance(x[1:2], x[3]))     # based on the ellipsoid
            create_backup()
            iBefehl                       <<- iBefehl+1
            befehl$updateClassification   <-  iBefehl
            befehl$updateVisualizations   <-  iBefehl
            befehl$updateEllipsoid        <-  iBefehl
            befehl$plot                   <-  iBefehl
          }
        }
      }
    })

    observe({
      tmpAxis1 = input$NumericMainAxis1
      disable("execEM")
      #print("Refresh Values for Main Axis 1")
      #print("AdaptDunes:Refresh Main Axis 1")
      if(!is.null(tmpAxis1)){
        if(is.numeric(tmpAxis1)){
          if(tmpAxis1 < MainAxesAngle[[CurrGauss]][2]){
            tmpVal = MainAxesAngle[[CurrGauss]][2]
            MainAxesAngle[[CurrGauss]][2] <<- tmpAxis1
            MainAxesAngle[[CurrGauss]][1] <<- tmpVal
            tmpVal = MainAxesAngle[[CurrGauss]][3]
            MainAxesAngle[[CurrGauss]][3] <<- MainAxesAngle[[CurrGauss]][4]
            MainAxesAngle[[CurrGauss]][4] <<- tmpVal
          }else{
            MainAxesAngle[[CurrGauss]][1] <<- tmpAxis1
          }
          updateNumericInput(session, "NumericMainAxis1", value = MainAxesAngle[[CurrGauss]][1])
          updateNumericInput(session, "NumericMainAxis2", value = MainAxesAngle[[CurrGauss]][2])
          updateNumericInput(session, "NumericAngle",     value = decimal2CompassDegree(MainAxesAngle[[CurrGauss]][3]))
          CovMats <<- lapply(MainAxesAngle,                                       # Update the Covariance Matrix
                             function(x) axesAngle2Covariance(x[1:2], x[3]))      # based on the ellipsoid
          create_backup()
          iBefehl                       <<- iBefehl+1
          befehl$updateNumericInputs    <-  iBefehl
          befehl$updateClassification   <-  iBefehl
          befehl$updateVisualizations   <-  iBefehl
          befehl$updateEllipsoid        <-  iBefehl
          befehl$plot                   <-  iBefehl
        }
      }
    })

    observe({
      tmpAxis2 = input$NumericMainAxis2
      disable("execEM")
      #print("Refresh Values for Main Axis 2")
      #print("AdaptDunes:Refresh Main Axis 2")
      if(!is.null(tmpAxis2)){
        if(is.numeric(tmpAxis2)){
          if(tmpAxis2 > MainAxesAngle[[CurrGauss]][1]){
            tmpVal = MainAxesAngle[[CurrGauss]][1]
            MainAxesAngle[[CurrGauss]][1] <<- tmpAxis2
            MainAxesAngle[[CurrGauss]][2] <<- tmpVal
            tmpVal = MainAxesAngle[[CurrGauss]][3]
            MainAxesAngle[[CurrGauss]][3] <<- MainAxesAngle[[CurrGauss]][4]
            MainAxesAngle[[CurrGauss]][4] <<- tmpVal
          }else{
            MainAxesAngle[[CurrGauss]][2] <<- tmpAxis2
          }
          updateNumericInput(session, "NumericMainAxis1", value = MainAxesAngle[[CurrGauss]][1])
          updateNumericInput(session, "NumericMainAxis2", value = MainAxesAngle[[CurrGauss]][2])
          updateNumericInput(session, "NumericAngle",     value = decimal2CompassDegree(MainAxesAngle[[CurrGauss]][3]))
          CovMats <<- lapply(MainAxesAngle,                                      # Update the Covariance Matrix
                             function(x) axesAngle2Covariance(x[1:2], x[3]))     # based on the ellipsoid
          create_backup()
          iBefehl                       <<- iBefehl+1
          befehl$updateNumericInputs    <-  iBefehl
          befehl$updateClassification   <-  iBefehl
          befehl$updateVisualizations   <-  iBefehl
          befehl$updateEllipsoid        <-  iBefehl
          befehl$plot                   <-  iBefehl
        }
      }
    })


    #--------------------------------------------------------------------------#
    #
    #
    #          Plotting
    #
    #
    #--------------------------------------------------------------------------#

    #--------------------------------------------------------------------------#
    # Plot view left
    #--------------------------------------------------------------------------#
    output$leftView <- renderPlotly({
      befehl$plot

      if(Debug){
        cat(file=stderr(), "Render Left Main Plotly was toggled","\n")
      }
      if(ControlPlotReaction1==F){
        ControlPlotReaction1<<-T
        if(values$leftPlot2 == "Model Classes"){
          plotClassification = OriginalCls[dataSampleIdx2]
        }else{
          plotClassification = Classification[dataSampleIdx2]
        }

        # leftPlot is button for switching to different state (so if Data Est. is shown, then PlotModelDotDensity3D is plotted)
        if(values$leftPlot == "Data Estimation"){
          CovMats <<- lapply(MainAxesAngle, function(x) axesAngle2Covariance(x[1:2], x[3]))
          if(values$DotOrGrid == "Grid Density"){
            plotOut = plotModelDotDensity3D(DataVis, XKernel, YKernel,
                                            EmpiricDataPDE[dataSampleIdx2], DotDensity,
                                            Means, CovMats, Weights, MainAxesAngle,
                                            gaussColors, plotClassification,
                                            cameraLeft,
                                            ShowScatter = values$dataScatter,
                                            ShowMarkers = values$ShowMarkersLeft,
                                            AxNames = AxNames,
                                            Source = "L1",
                                            Debug = Debug)
          }else{
            plotOut = plotModelGridDensity3D(DataVis,
                                             XKernel, YKernel,
                                             ContinuousDataPDE,
                                             GridDensity,
                                             Means, CovMats, Weights,
                                             MainAxesAngle,
                                             gaussColors, plotClassification,
                                             cameraLeft,
                                             AxNames = AxNames,
                                             Source = "L1",
                                             Debug = Debug)
          }
        }else if(values$leftPlot == "Model PDF"){
          plotOut = plotDataGridDensity3D(DataVis, XKernel, YKernel,
                                          ContinuousDataPDE, EmpiricDataPDE[dataSampleIdx2],
                                          #Means, CovMats, MainAxesAngle,
                                          gaussColors, plotClassification, cameraRight,
                                          ShowScatter = values$dataScatter,
                                          AxNames = AxNames,
                                          Source = "L1",
                                          Debug = F)
        }
        plotOut
      }
    })

    #--------------------------------------------------------------------------#
    # Plot view right
    #--------------------------------------------------------------------------#
    output$rightView  = renderPlotly({
      befehl$plot

      if(Debug){
        cat(file=stderr(), "Render Right Main Plotly:","\n")
      }

      if(ControlPlotReaction2==F){
        ControlPlotReaction2<<-T
        if(values$showAxis == "Axis"){
          ShowAxis = FALSE
        }else{
          ShowAxis = TRUE
        }
        if(values$showEllipsoids == "Ellipsoids"){
          ShowEllipsoids = FALSE
        }else{
          ShowEllipsoids = TRUE
        }
        if(values$showScatter == "Scatter"){
          ShowScatter = FALSE
        }else{
          ShowScatter = TRUE
        }
        if(values$rightCls == "Model Classes"){
          plotClassification = OriginalCls[dataSampleIdx2]
        }else{
          plotClassification = Classification[dataSampleIdx2]
        }
        shapes       = isolate(values$shapes)
        shapeText    = isolate(values$shapeText)
        CovMats <<- lapply(MainAxesAngle, function(x) axesAngle2Covariance(x[1:2], x[3]))
        if(values$rightPlot == 1){
          plotOut = plotScatter(Data = DataVis, CurrGauss = CurrGauss,
                                Means = Means, Covariances = CovMats, MainAxesAngle = MainAxesAngle,
                                Colors = gaussColors, Cls = plotClassification,
                                Shapes = shapes, ShapeText = shapeText, AxNames = AxNames,
                                ShowAxis = ShowAxis, ShowEllipsoids = ShowEllipsoids,
                                ShowGaussNr = FALSE, Source = "D")
        }else if(values$rightPlot == 2){
          plotOut = plotDensity(Data = DataVis, Cls = plotClassification, CurrGauss = CurrGauss, Colors = gaussColors,
                                GridDensity = GridDensity,
                                Means = Means, Covariances = CovMats, Weights = Weights, MainAxesAngle = MainAxesAngle,
                                XKernel = XKernel, YKernel = YKernel,
                                Shapes = shapes, ShapeText = shapeText, AxNames = AxNames,
                                ShowAxis = ShowAxis, ShowEllipsoids = ShowEllipsoids,
                                ShowScatter = ShowScatter,
                                ShowGaussNr = FALSE, Source = "D")
        }else if(values$rightPlot == 3){
          plotOut = plotDataDensity(Data = DataVis, Cls = plotClassification, CurrGauss = CurrGauss, Colors = gaussColors,
                                    ContinuousDataPDE = ContinuousDataPDE,
                                    Means = Means, Covariances = CovMats, MainAxesAngle = MainAxesAngle,
                                    XKernel = XKernel, YKernel = YKernel,
                                    Shapes = shapes, ShapeText = shapeText, AxNames = AxNames,
                                    ShowAxis = ShowAxis, ShowEllipsoids = ShowEllipsoids,
                                    ShowScatter = ShowScatter,
                                    ShowGaussNr = FALSE, Source = "D")
        }else{
          plotOut = plotPoliticalMap(Data = DataVis, CurrGauss = CurrGauss,
                                     Means = Means, Covariances = CovMats, Weights = Weights, MainAxesAngle = MainAxesAngle,
                                     Colors = gaussColors, Cls = plotClassification,
                                     Shapes = shapes, ShapeText = shapeText, AxNames = AxNames,
                                     ShowAxis = ShowAxis, ShowEllipsoids = ShowEllipsoids,
                                     ShowScatter = ShowScatter, Source = "D")
        }
        plotOut
      }
    })

    #--------------------------------------------------------------------------#
    #                   Update and Manage Ellipsoids
    #--------------------------------------------------------------------------#
    observe({
      befehl$updateEllipsoid
      if(Debug){
        cat(file = stderr(), "Update Ellipsoids\n")
      }
      shapes     = list()
      shapeText  = matrix(ncol = 3, nrow = 0)
      Annotation = rbind()
      tryCatch({
        CovMats <<- lapply(MainAxesAngle, function(x) axesAngle2Covariance(x[1:2], x[3]))
        for(i in 1:length(Means)){
          el        = mixtools::ellipse(Means[[i]], CovMats[[i]], draw = F)
          svgCPrep  = paste0("", paste0(apply(el, 1, function(x) paste0("L ",paste0(round(x,4), collapse=" "))),
                                        collapse = " "))
          svgC      = paste0("M",substr(svgCPrep, 2, nchar(svgCPrep)))
          ellipse   = list(type = 'path', fillcolor = gaussColors[i], opacity = 0.3, path = svgC)
          shapes    = c(shapes, list(ellipse))#, vec1, vec2))
          shapeText = rbind(shapeText, c(Means[[i]][1], Means[[i]][2], i))
          values$lastKnownShapeStart[i,] = round(el[i,],4)

          values$shapes    = shapes
          values$shapeText = shapeText
        }
      },error =  function(e){cat(file = stderr(), "Some  error happened:"); cat(file = stderr(), paste0(e))})
    })

    #--------------------------------------------------------------------------#
    #                  React to shape movements
    #--------------------------------------------------------------------------#
    #   reacts: event_data("plotly_relayout")
    #   triggers:  numericInputs,   values$gExp
    observe({
      if(Debug){
        cat(file=stderr(), "Plotly Shape Event is triggered:", "\n")
      }
      ed = event_data("plotly_relayout", source = "D")    # Only reacts on events from plots of source "D"
      #browser()
      if(!is.null(ed) & !is.null(Means)){
        x = names(ed)[1]
        #browser()
        preR = gregexpr("shapes\\[(\\d+)\\].path",x, perl = T)
        if(preR[[1]] != -1){
          if(Debug) cat(file=stderr(), "Plotly Event contains a Shape change exp from: ")
          rmatches  = preR[[1]]
          capSt     = attr(rmatches, 'capture.start')
          capLength = attr(rmatches, 'capture.length')
          captures  = substring(x, c(capSt), c(capSt)+c(capLength)-1)
          shapeNr   = as.integer(captures) + 1
          if(Debug) cat(file=stderr(), as.character(Means[[shapeNr]]))
          if(Debug) cat(file=stderr(), " to ")
          #         browser()
          a      = strsplit(ed[[1]], " ")[[1]][2:3]
          newPos = round(as.numeric(a),4)
          mov    = isolate(values$lastKnownShapeStart[shapeNr,]) - newPos # Check if movement happened
          if(all(mov == c(0,0))){                                         # Difference of positions is 0 in case of no movements
            cat(file=stderr(), " (no actual change happened) ")
          }
          else{
            #print(paste0("Updating Gaussian."))
            Means[[shapeNr]] <<- Means[[shapeNr]] - mov
            create_backup()
            #values$updateEllipsoid = values$updateEllipsoid + 1
            iBefehl                       <<- iBefehl+1
            befehl$updateVisualizations   <-  iBefehl
            befehl$updateEllipsoid        <-  iBefehl
            befehl$updateNumericInputs    <-  iBefehl
            #befehl$updateClassification   <-  iBefehl
            #befehl$updateSlider           <-  iBefehl
            #befehl$plot                   <-  iBefehl
          }
          if(Debug) cat(file=stderr(), as.character(Means[[shapeNr]]))
          if(Debug) cat(file=stderr(), " for ")
          if(Debug) cat(file=stderr(), as.character(shapeNr))
          if(Debug) cat(file=stderr(), "\n")
        }
      }
    })

    #--------------------------------------------------------------------------#
    #  Remember camera angle of 3D plots
    #--------------------------------------------------------------------------#
    # Capture camera angle of 3D plots and keep it through changes of other
    # paramters and update camera angle on change by user
    # Track changes via plotly events
    observe({
      if(Debug){
        cat(file=stderr(), "Plotly 3D Event is triggered:", "\n")
      }
      ed = event_data("plotly_relayout", source = "L1")
      if(!is.null(ed)){
        cameraLeft <<- ed$scene.camera
      }
    })

    observe({
      if(Debug){
        cat(file=stderr(), "Plotly 3D Event is triggered:", "\n")
      }
      ed = event_data("plotly_relayout", source = "L2")
      if(!is.null(ed)){
        cameraLeft <<- ed$scene.camera
      }
    })

    observe({
      ed = event_data("plotly_relayout", source = "R")
      if(!is.null(ed)){
        cameraRight <<- ed$scene.camera
      }
    })

    #--------------------------------------------------------------------------#
    #  Select Gaussian Component
    #--------------------------------------------------------------------------#
    # Change current Gauss
    observe({
      #print("AdaptDunes: Update CurrGauss")
      #print("Update CurrGauss")
      if (!is.null(input$gaussianNo)){
        CurrGauss                     <<- input$gaussianNo
        iBefehl                       <<- iBefehl+1
        befehl$updateCurrGauss        <-  iBefehl
        befehl$drawGaussianComponents <-  iBefehl
        befehl$drawUIControl <- iBefehl
      }
    })

    #--------------------------------------------------------------------------#
    #  Load classes and save settings
    #--------------------------------------------------------------------------#
    observeEvent(input$loadCls,{
      tmpFile = rstudioapi::selectFile()
      if(!is.null(tmpFile)){
        # Read classification (either 'cls' or 'csv' format)
        if(substr(tmpFile, nchar(tmpFile) - 2, nchar(tmpFile)) == "cls"){
          tmpV   = dbt.DataIO::ReadCLS(FileName = tmpFile)
          newCls = tmpV$Cls
        }else if (substr(tmpFile, nchar(tmpFile) - 2, nchar(tmpFile)) == "csv"){
          tmpV   = read.delim2(file = tmpFile, header = FALSE)
          if(is.list(tmpV)){
            tmpV = unlist(tmpV)
          }
          if(is.na(as.numeric(tmpV[1]))){
            newCls = as.numeric(tmpV[2:length(tmpV)])
          }else{
            newCls = as.numeric(tmpV)
          }
        }else{
          newCls = 0
        }
        # Check if exact number of entries are there
        if(is.numeric(newCls) & is.vector(newCls) & length(newCls) == length(OriginalCls)){
          OriginalCls <<- newCls
          iBefehl              <<- iBefehl+1
          befehl$buttonClasses <- iBefehl
        }else{
          cat(file=stderr(), "Given class has not the correct number of class labels.", "\n")
        }
      }
    })

    observeEvent(input$saveAdaptDunes,{
      tmpDirBackUp = getwd()
      tmpDir = rstudioapi::selectDirectory()
      if(!is.null(tmpDir)){
        setwd(tmpDir)
        tmpTime = gsub(":", "", format(Sys.time(), "%Y%m%d%X"))

        tmpSamDat = Data[, c(FeatureIdx1, FeatureIdx2)]
        tmpV      = computeGMMClassification(tmpSamDat, Means, CovMats, Weights)
        Classification <<- tmpV$Classification

        #saveRDS(object = OriginalCls, file = paste0(tmpTime, "_Original_Classes.Rds"))
        #saveRDS(object = Classification, file = paste0(tmpTime, "_Model_Classes.Rds"))
        SaveMeans   = round(matrix(unlist(Means), byrow = TRUE, ncol = 2), 2)
        SaveCovMats = round(matrix(unlist(lapply(CovMats, t)), byrow = TRUE, ncol = 2), 2)
        SaveWeights = round(Weights, 2)
        SaveMAA     = round(matrix(unlist(MainAxesAngle), byrow = TRUE, ncol = 4), 2)

        if(dbt == TRUE){
          dbt.DataIO::WriteLRN(FileName = "Means",         Data = SaveMeans,   Key = 1:length(Means), Header = AxNames, CommentOrDigits = "Mean 1x2 stacked vertically.")
          dbt.DataIO::WriteLRN(FileName = "CovMats",       Data = SaveCovMats, Key = 1:(length(Means)*2), Header = AxNames, CommentOrDigits = "Covariance Matrix 2x2 stacked vertically.")
          dbt.DataIO::WriteLRN(FileName = "Weights",       Data = SaveWeights, Key = 1:length(Means), Header = "GMMWeights")
          dbt.DataIO::WriteLRN(FileName = "MainAxesAngle", Data = SaveMAA,     Key = 1:length(Means), Header = c("Ax1", "Ax2", "AngleForAx1", "AngleForAx2"), CommentOrDigits = "Axes information for all GMM components stacked vertically.")
          dbt.DataIO::WriteCLS(FileName = "ModelClassification", Cls = Classification, Key = 1:length(Cls))
        }else{
          write.table(x = SaveMeans,      file = "Means.csv",               row.names = FALSE, col.names = AxNames)
          write.table(x = SaveCovMats,    file = "CovMats.csv",             row.names = FALSE, col.names = AxNames)
          write.table(x = SaveWeights,    file = "Weights.csv",             row.names = FALSE, col.names = "GMMWeights")
          write.table(x = SaveMAA,        file = "MainAxesAngle.csv",       row.names = FALSE, col.names = c("Ax1", "Ax2", "AngleForAx1", "AngleForAx2"))
          write.table(x = Classification, file = "ModelClassification.csv", row.names = FALSE, col.names = FALSE)
        }
        #save(list = c("Means", "CovMats", "Weights", "MainAxesAngle",
        #              "Classification", "OriginalCls"),
        #     file = paste0(tmpTime, "_AdaptDunes.Rdata"))
        setwd(tmpDirBackUp)
      }
    })

    #--------------------------------------------------------------------------#
    #  Add data scatter or markers
    #--------------------------------------------------------------------------#
    observeEvent(input$dataScatter,{
      values$dataScatter = !values$dataScatter
    })

    observeEvent(input$showMarkersLeft,{
      values$ShowMarkersLeft = input$showMarkersLeft
    })

    #--------------------------------------------------------------------------#
    #  Change Plot of left side (Model PDF or Data Estimation)
    #--------------------------------------------------------------------------#
    observeEvent(input$ChangeDotOrGrid,{
      if(values$DotOrGrid == "Dot Density"){
        if(Debug){
          print("Change to the Grid Density")
        }
        values$DotOrGrid = "Grid Density"
      }else{
        if(Debug){
          print("Change to the Dot Density")
        }
        values$DotOrGrid = "Dot Density"
      }
    })

    observeEvent(input$changeLeftPlot,{
      if(values$leftPlot == "Data Estimation"){
        if(Debug){
          print("Change to the Model PDF.")
        }
        values$leftPlot = "Model PDF"
      }else{
        if(Debug){
          print("Change to the Date Estimation.")
        }
        values$leftPlot = "Data Estimation"
      }
    })

    observeEvent(input$changeLeftPlot2,{
      if(values$leftPlot2 == "Orig. Classes"){
        if(Debug){
          print("Change to the Model Classes.")
        }
        values$leftPlot2 = "Model Classes"
      }else{
        if(Debug){
          print("Change to the Original Classes.")
        }
        values$leftPlot2 = "Orig. Classes"
      }
    })

    #--------------------------------------------------------------------------#
    #  Change Plot of right side (Scatter Plot or Classification Map)
    #--------------------------------------------------------------------------#
    # Upper buttons for control of right plots
    observeEvent(input$changeRightPlot1,{
      if(RightPlotControlStopReactions==F){
        RightPlotControlStopReactions<<-T
        values$rightPlot = 1
      }
    })

    observeEvent(input$changeRightPlot2,{
      if(RightPlotControlStopReactions==F){
        RightPlotControlStopReactions<<-T
        values$rightPlot = 2
      }
    })

    observeEvent(input$changeRightPlot3,{
      if(RightPlotControlStopReactions==F){
        RightPlotControlStopReactions<<-T
        values$rightPlot = 3
      }
    })

    observeEvent(input$changeRightPlot4,{
      if(RightPlotControlStopReactions==F){
        RightPlotControlStopReactions<<-T
        values$rightPlot = 4
      }
    })

    # Lower buttons for control of right plots
    observeEvent(input$showAxis,{
      if(RightPlotControlStopReactions==F){
        if(values$showAxis == "No axis"){
          if(Debug){
            print("Showing no Axis")
          }
          values$showAxis = "Axis"
        }else{
          if(Debug){
            print("Axis")
          }
          values$showAxis = "No axis"
        }
      }
    })

    observeEvent(input$showEllipsoids,{
      if(RightPlotControlStopReactions==F){
        if(values$showEllipsoids == "No ellipsoids"){
          if(Debug){
            print("Showing no ellipsoids")
          }
          values$showEllipsoids = "Ellipsoids"
        }else{
          if(Debug){
            print("Ellipsoids")
          }
          values$showEllipsoids = "No ellipsoids"
        }
      }
    })

    observeEvent(input$showScatter,{
      if(RightPlotControlStopReactions==F){
        if(values$showScatter == "No scatter"){
          if(Debug){
            print("Showing no scatter")
          }
          values$showScatter = "Scatter"
        }else{
          if(Debug){
            print("Showing scatter")
          }
          values$showScatter = "No scatter"
        }
      }
    })

    observeEvent(input$changeRightCls,{
      if(RightPlotControlStopReactions==F){
        if(values$rightCls == "Orig. Classes"){
          values$rightCls = "Model Classes"
        }else{
          values$rightCls = "Orig. Classes"
        }
      }
    })

    #--------------------------------------------------------------------------#
    #  Add/Remove Gaussian Components
    #--------------------------------------------------------------------------#
    observeEvent(input$addGaussButton,{
      if(ControlAddGauss==F){
        ControlAddGauss<<-T
        if(Debug){
          cat(file=stderr(), "Added Gauss:", "\n")
        }
        NumGauss                  <<- NumGauss + 1
        Means[[NumGauss]]         <<- stats::runif(2, -1, 1)
        CovMats[[NumGauss]]       <<- diag(2)
        MainAxesAngle[[NumGauss]] <<- c(1,1,0,90) # Unit length for axes (1,1), 0 angle of 1st axis 1, 90 angle of 2nd axis
        Weights                   <<- rep(1/NumGauss, NumGauss)
        CurrGauss                 <<- NumGauss

        iBefehl                       <<- iBefehl+1
        befehl$updateClassification   <-  iBefehl
        befehl$updateVisualizations   <-  iBefehl
        befehl$updateCurrGauss        <-  iBefehl
        befehl$drawUIControl          <-  iBefehl
        befehl$drawGaussianComponents <-  iBefehl
        befehl$updateGaussianNo       <-  iBefehl
        befehl$updateNumericInputs    <-  iBefehl
        #befehl$updateSlider           <-  iBefehl
        #befehl$plot                   <-  iBefehl
        befehl$updateEllipsoid        <-  iBefehl

        if(Debug) cat(file=stderr(), as.character(Weights), "\n")
      }
    })

    observeEvent(input$remGaussButton,{
      if(ControlRemoveGauss==F){
        ControlRemoveGauss<<-T
        if(length(Means) > 1){
          print("Removing Gaussian.")
          if(CurrGauss < NumGauss){                          # If current shown Gaussian component is less
            TmpColor = gaussColors[CurrGauss]                # than the number of components, then move
            for(i in CurrGauss:(NumGauss-1)){                # the components in the list one forward in
              Means[[i]]         <<- Means[[i+1]]            # order to remove the current one
              CovMats[[i]]       <<- CovMats[[i+1]]
              MainAxesAngle[[i]] <<- MainAxesAngle[[i+1]]
              Weights[i]         <<- Weights[i+1]
              gaussColors[i]     <<- gaussColors[i+1]
            }
            gaussColors[NumGauss] <<- TmpColor
          }
          if(CurrGauss==NumGauss){                           # If the current shown Gaussian (to be removed) is
            CurrGauss <<- NumGauss-1                         # the last component, move one forward
          }
          NumGauss      <<- NumGauss - 1                     # The number of Gaussians is decimated by one
          Means         <<- Means[1:NumGauss]                # Use the new number of Gaussian to assign the new list sizes
          CovMats       <<- CovMats[1:NumGauss]
          MainAxesAngle <<- MainAxesAngle[1:NumGauss]
          Weights       <<- rep(1 / NumGauss, NumGauss)

          iBefehl                       <<- iBefehl+1        # Trigger next command
          befehl$updateClassification   <-  iBefehl          # Start with internal changes such as classification first
          befehl$updateVisualizations   <-  iBefehl
          befehl$updateCurrGauss        <-  iBefehl          # Continue drawings
          #befehl$updateSlider           <-  iBefehl
          befehl$updateNumericInputs    <-  iBefehl
          befehl$updateEllipsoid        <-  iBefehl
          befehl$drawUIControl          <-  iBefehl
          befehl$drawGaussianComponents <-  iBefehl
          befehl$updateGaussianNo       <-  iBefehl
          #befehl$plot                   <-  iBefehl          # Toggel plots last
        }else{
          cat(file = stderr(), "Last Gaussian component cannot be removed. \n")
        }
      }
    })

    #--------------------------------------------------------------------------#
    #         Expectation-Maximization Algorithm Interaction
    #--------------------------------------------------------------------------#
    observeEvent(input$execEM,{
      disable("execEM")
      if(ControlEM==F){
        ControlEM<<-T
        print("Start EM")
        #browser()
        tryCatch({
          #if(length(Weights) < 2){
          #  cat(file = stderr(), "There must be at least 2 Gaussian mixture components.\n")
          #  return()
          #}
          create_backup()
          res = ad_interfaceEM(Data     = Data[ ,c(FeatureIdx1, FeatureIdx2)],
                               Mu       = Means,
                               Sigma    = CovMats,
                               Lambda   = Weights,
                               Addmodes = EMSplitModes,
                               Nit      = EMNumIterations)
          Means         <<- res$Mu
          CovMats       <<- res$Sigma
          Weights       <<- res$Lambda
          MainAxesAngle <<- res$MainAxesAngle
          NumGauss      <<- length(Means)
          if(Debug){
            cat(file = stderr(), "EM finished \n")
          }
          iBefehl                       <<- iBefehl+1
          befehl$updateVisualizations   <-  iBefehl
          #befehl$updateSlider           <-  iBefehl
          befehl$updateNumericInputs    <-  iBefehl
          befehl$updateCurrGauss        <-  iBefehl
          befehl$updateEllipsoid        <-  iBefehl
          befehl$drawUIControl          <-  iBefehl
          befehl$drawGaussianComponents <-  iBefehl
          #befehl$plot                   <-  iBefehl
        },
         error = function(e) cat(file = stderr(), paste0("EM: Error", e, "\n")),
         warning = function(e) cat(file = stderr(), paste0("EM: Error", e, "\n"))
        )
      }
    })

    observeEvent(input$undoEM,{
      if(ControlUndoRedoEM==F){
        ControlUndoRedoEM<<-T
        print("Undo")
        # browser()
        recover_latest_backup()
      }
    })

    observeEvent(input$redoEM,{
      if(ControlUndoRedoEM==F){
        ControlUndoRedoEM<<-T
        print("Redo")
        # browser()
        redo_next_backup()
      }
    })

    #--------------------------------------------------------------------------#
    #  Norm Weights (All or Others)
    #--------------------------------------------------------------------------#
    observeEvent(input$normAll,{
      if(ControlNormAll==F){
        ControlNormAll<<-T
        SumWeights =   sum(Weights)
        Weights    <<- Weights / SumWeights
        #cat(paste0("Norm All:", Weights))
        if(Debug){
          cat(file = stderr(), paste0("Norm all GMM components: ", Weights, "\n"))
        }
        iBefehl                <<- iBefehl+1
        #befehl$updateSlider    <-  iBefehl
        befehl$updateNumericInputs    <-  iBefehl
        #befehl$plot            <-  iBefehl
        befehl$updateVisualizations   <-  iBefehl
        befehl$updateEllipsoid <-  iBefehl
      }
    })

    observeEvent(input$normOth,{
      if(ControlNormOthers==F){
        ControlNormOthers<<-T
        SumWeight       = sum(Weights)
        SumWeightOthers = SumWeight - Weights[CurrGauss]
        TargetSize      = 1 - Weights[CurrGauss]
        for(i in 1:NumGauss){
          if(i!=CurrGauss){
            Weights[i] <<- Weights[i]/SumWeightOthers*TargetSize # SumWeightOthers/(NumGauss-1)
          }
        }
        if(Debug){
          cat(file = stderr(), paste0("Norm other GMM components: ", Weights, "\n"))
        }
        iBefehl                       <<- iBefehl+1
        #befehl$updateSlider           <-  iBefehl
        befehl$updateNumericInputs    <-  iBefehl
        #befehl$plot                   <-  iBefehl
        befehl$updateVisualizations   <-  iBefehl
        befehl$updateEllipsoid        <-  iBefehl
      }
    })

    #--------------------------------------------------------------------------#
    #               End: Quit, Close
    #--------------------------------------------------------------------------#
    observeEvent(input$quitButton,{
      print("Close App")
      V                      = computeGMMClassification(Data[,c(FeatureIdx1, FeatureIdx2)],
                                                        Means, CovMats, Weights)
      Classification         <<- V$Classification
      stopApp(list("Means"                   =   Means,
                   "CovarianceMatrices"      =   CovMats,
                   "Weights"                 =   Weights,
                   "PrincipalComponentAxis"  =   MainAxesAngle[1:2],
                   "Angle"                   =   MainAxesAngle[3],
                   "Cls"                     =   Classification))
    })

    session$onSessionEnded(function() {
      print("Close App")
      V                      = computeGMMClassification(Data[,c(FeatureIdx1, FeatureIdx2)],
                                                        Means, CovMats, Weights)
      Classification         <<- V$Classification
      stopApp(list("Means"                   =   Means,
                   "CovarianceMatrices"      =   CovMats,
                   "Weights"                 =   Weights,
                   "PrincipalComponentAxis"  =   MainAxesAngle[1:2],
                   "Angle"                   =   MainAxesAngle[3],
                   "Cls"                     =   Classification))
    })
  }
  return(shiny::runApp(shinyApp(ui,server)))
}
