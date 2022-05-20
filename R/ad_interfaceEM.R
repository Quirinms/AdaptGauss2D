ad_interfaceEM = function(Data, Mu, Sigma, Lambda, Nmodes = 1, Addmodes = 1,
                          Nit = 1, Verbose = 0){
  # DESCRIPTION
  # Interface to Expectation Maximization algorithm from Dr. Baggenstoss.
  #
  # INPUT
  # Data[1:n, 1:2]    Numeric matrix with n observations and 2 features.
  # Mu                List with l [1:2] numerical vector defining the means of
  #                   the l GMM components.
  # Sigma             List with l [1:2, 1:2] numerical matrices defining the
  #                   covariance matrices of the l GMM components.
  # Lambda[1:l]       Numerical vector with weights for each GMM component.
  # Nit               Integer stating the number of iterations.
  # Addmodes          0 or 1. 1: Allow EM to add modes. 0: forbid adding modes.
  # Verbose           Integer indicating how much information output should be given (0,1,2).
  #
  # OUTPUT
  # Res    List of GMM settings computed by an EM algorithm of choice.
  #
  # DETAILS
  # If you use EM Baggenstoss for initialization, use many iterations (100-500)
  # If you use EM Baggenstoss to optimize manual adjustments, try 1 iteration
  #
  # Author: QMS 15.12.2021

  if(!is.matrix(Data)){
    print("Data is not a matrix.")
    return()
  }
  if(!is.numeric(Data)){
    print("Data is not a numeric matrix.")
    return()
  }

  if(!is.numeric(Nit)){
    Nit = 1
  }

  if(is.null(Mu) | is.null(Sigma) |is.null(Lambda)){
    Mu     = NULL
    Sigma  = NULL
    Lambda = NULL
  }

  if(!is.null(Mu)){
    if(!is.list(Mu)){
      print("Mu is not a list.")
      return()
    }else{
      for(i in 1:length(Mu)){
        if(!is.vector(Mu[[i]]) & !is.numeric(Mu[[i]])){
          print("Mu is not a list of numeric vectors.")
          return()
        }
      }
    }

    if(!is.list(Sigma)){
      print("Sigma is not a list.")
      return()
    }else{
      for(i in 1:length(Sigma)){
        if(!is.matrix(Sigma[[i]]) & !is.numeric(Sigma[[i]])){
          print("Sigma is not a list of numeric matrices.")
          return()
        }
      }
    }

    if(!is.vector(Lambda)){
      print("Lambda is not a vector.")
      return()
    }
    if(!is.numeric(Lambda)){
      print("Lambda is not a numeric matrix.")
      return()
    }
  }
  MinStd = rep(0.1, dim(Data)[2])
  if(!is.null(Mu) & !is.null(Sigma) & !is.null(Lambda)){
    nrGauss    = length(Mu)
    print("Executing Baggenstoss EM Algorithm")
    Names      = unlist(lapply(1:dim(Data)[2], function(x) paste0("GMM",x)))
    Covariance = do.call(rbind, Sigma)
    Mean       = do.call(rbind, Mu)
    Weight     = Lambda
    Parm       = list("Features" = list("Names"  = Names,
                                        "MinStd" = MinStd),
                      "Modes"  = list("Covariance" = Covariance,
                                      "Mean"       = Mean,
                                      "Weight"     = Lambda))
  }else{
    print("Initialization with Baggenstoss EM Algorithm")
    InitEM  = init_gmix(Data, Nmodes = Nmodes, MinStd = MinStd, Names=NULL, RandInit=1, Verbose=0)
    Parm    = InitEM$Parm
  }
  TrainEM = ad_gmix_trainscript(Parm = Parm, Data = Data, Nit = Nit,
                                Addmodes = Addmodes, Verbose = Verbose)
  Mu      = lapply(seq_len(nrow(TrainEM$EMmean)), function(i) TrainEM$EMmean[i,])
  Lambda  = as.vector(TrainEM$EMalpha)
  CovMat  = lapply(1:(nrow(TrainEM$EMCov)/2), function(x) TrainEM$EMCov[(2*(x-1)+1):(2*x),])
  # Covariance matrix to principal components axes lengths and angle
  MainAxesAngle = lapply(CovMat, function(x) covariance2AxesAngle(x))
  # Baggenstoss: Expectation Maximization from Dr. Baggenstoss does not always
  # yield symmetric postive definite covariance matrices
  # Quirin => Engineered solution: use SVD + ellipsoidal computation to ensure
  # positive definite covariance matrices at all times
  CovMat = lapply(MainAxesAngle, function(x) axesAngle2Covariance(x[1:2], x[3]))
  MainAxesAngle = lapply(CovMat, function(x) covariance2AxesAngle(x))
  #MainAxesAngle = lapply(MainAxesAngle, function(x) KeepCirclesElliptic(x))
  #CovMat = lapply(MainAxesAngle, function(x) axesAngle2Covariance(x[1:2], x[3]))
  return(list("Mu"            = Mu,
              "Sigma"         = CovMat,
              "Lambda"        = Lambda,
              "MainAxesAngle" = MainAxesAngle))
}
