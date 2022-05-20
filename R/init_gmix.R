init_gmix = function(Data, Nmodes, MinStd, Names=NULL, RandInit=1, Verbose=0){
  # V = init_gmix(Data, Nmodes, MinStd, Names, RandInit=1, Verbose=1)
  # V$Parm
  # 
  # DESCRIPTION
  # Initialization of parameters for Baggenstoss EM alg.
  # 
  # INPUT
  # Data[1:n, 1:d]     Numerical matrix with n samples and d feature dimensions.
  #                    Matlab implemntation commentary: normalized data.
  # Nmodes[1:l]        Number of modes to initialize.
  # MinStd[1:d]        Numerical vector with covariance constraints for each
  #                    feature. Preventing the covariance matrix to become
  #                    singular.
  # Names[1:d]         String vector with d feature Names
  # 
  # OPTIONAL
  # RandInit           Use randomized initialization. Random=1, Nonrandom=0.
  # Verbose            Verbose = 0: do not print messages
  #                    Verbose = 1: print messages
  #                    Default: Verbose = 1
  # 
  # OUTPUT
  # Parm             Nested list with parameters for GMM.
  #                  Features carrying permanent values:
  #                  Features$Names [1:d]      String vector with feature names.
  #                  Features$MinStd [1:l]     Vector of covariance constraints.
  #                  Modes list of modifyable GMM parameters for all modes.
  #                  Modes$Covariance [1:d*l, 1:d]     Numerical matrix with l
  #                        many square matrices stacked vertically with the
  #                        covariance matrix.
  #                  Modes$Mean [1:l, 1:d]             Numerical matrix with l
  #                        row vectors containing the d dimensional means.
  #                  Modes$Weight [1:l]                Numerical vector Modes
  #                        weights for each mean.
  # 
  # DETAILS
  # The data has n samples and d Features There are l modes.
  # For each mode there is a weight. Each mode is d dimensional.
  # Each mode has one mean, which is a d dimensional vector.
  # Covariance denotes the covariance matrices of all modes.
  # See Matlab reference for more information.
  # 
  # 
  # 
  # in /dbt/EMforGauss/
  # Author QS 2021
  #----------------------------------------------------------------------------#
  # Check if required input variables are there with correct data structure
  if(is.null(Data)){
    stop("Input parameter Data is not given. Returning.")
  }
  if(!is.matrix(Data)){
    stop("Data is not a matrix!")
  }
  if(is.null(Nmodes)){
    stop("Provide the number of modes to start with!")
  }
  if(is.na(as.integer(Nmodes))){
    stop("Nmodes must be an integer!")
  }
  Nmodes = as.integer(Nmodes)
  if(is.null(MinStd)){
    stop("Provide covariance constraints!")
  }
  if(is.null(MinStd)){
    stop("Provide covariance constraints!")
  }
  if(!is.vector(MinStd)){
    stop("MinStd must be a numerical vector!")
  }
  if(!is.numeric(MinStd)){
    stop("MinStd must be a numerical vector!")
  }
  #----------------------------------------------------------------------------#
  # Check if there are any Names to work with, since there must be Names
  # either as input or as colnames in Data!
  if(is.null(Names)){
    if(is.null(colnames(Data))){
      stop("Names and colnames(Data) is NULL. Provide at least one of them!")
    }
    Names = colnames(Data)
  }
  #----------------------------------------------------------------------------#
  n_samples = dim(Data)[1]
  n_feats   = dim(Data)[2]
  Parm = list()
  Parm$Features$Names    = Names
  Parm$Features$MinStd  = MinStd
  # Determine mean and std of input Data
  if(n_samples > 1){
    data_means = apply(Data,2,function(x) mean(na.omit(x)))
    data_std   = apply(Data,2,function(x) sd(na.omit(x))) # Normalize the result by N-1 (built-in stats::sd functionality)
  }else{
    data_means = Data
    data_std   = matrix(0, n_feats, 1)
  }
  if(any((abs(data_means) > 1000*data_std)) || (max(data_std) > 100*min(data_std))){
    warning("init_gmix GSUM_INIT: Red Alert!! Warning: Possible ill-conditioning:")
    warning("init_gmix: Your features are shit. You should scale them or remove means!")
  }
  starting_std = 0 # Init
  if(n_samples > 1){
    xmin = apply(Data,2,function(x) min(na.omit(x)))
    xmax = apply(Data,2,function(x) max(na.omit(x)))
  }else{
    xmin = Data
    xmax = Data
  }
  # The initial value of covariances
  #starting_std = nanmax((xmax-xmin)/4,nanmax(MinStd(:)));
  starting_std = pmax((xmax-xmin)/4, max(c(NA, NA)), na.rm = TRUE)
  ##### ALU changed form: starting_std = max((xmax-xmin)/4,MinStd(:));
  # ...... see matlab file!
  if(RandInit){
    #--- Select starting means from a random set of data
    #--- Make sure there are none the same
    idx_has_duplicates = 1
    n_modes_min = 1
    while(idx_has_duplicates & Nmodes >= n_modes_min){
      i_tries = 0
      # Try 1000 times
      while(idx_has_duplicates & i_tries < 1000){
        # Generate Nmodes random numbers in [1,n_samples]
        idx = 1 + floor(pracma::rand(1,Nmodes) * n_samples)
        idx = sort(idx)
        # See if they are each unique
        if(length(idx) > 1){
          if(min(idx[2:Nmodes] - idx[1:Nmodes-1]) > 0){
            idx_has_duplicates = 0
          }
        }else{
          idx_has_duplicates = 0
        }
        i_tries = i_tries + 1
      }
      # if this still fails after 1000 tries, reduce number of modes and keep trying
      if(idx_has_duplicates){
        Nmodes = Nmodes - 1
        warning('gmix_init: reducing Nmodes to ', Nmodes)
        warning('gmix_init:  no. of samples ', n_samples, ' too small compared to Nmodes')
      }
    }
    # Number of modes has fallen below minimum, fail
    if(idx_has_duplicates){
      stop('gmix_init: Nmodes has been reduced to less than allowed, you need more data')
    }
  }else{ # non-random mode mean initialization
    Nmodes = min(Nmodes,n_samples)
    idx = 1:Nmodes
  }
  if(length(idx)==1){
    means = t(as.matrix(Data[idx,]))
  }else{
    means = Data[idx,]
  }
  # Create the (Cholesky of) covariances and store the means
  wts = matrix(1, 1, Nmodes)/Nmodes
  
  cholesky_covar = rbind()
  mean           = rbind()
  weight         = rbind()
  for(i_mode in 1:Nmodes){
    dim1 = dim(as.matrix(starting_std))[1]
    # Save results combined in one data structure
    cholesky_covar = rbind(cholesky_covar, diag(starting_std, dim1, dim1))
    mean           = rbind(mean, means[i_mode,])
    weight         = rbind(weight, wts[i_mode])
  }
  Parm$Modes$Covariance     = cholesky_covar
  Parm$Modes$Mean           = mean
  Parm$Modes$Weight         = weight
  if(Verbose > 0){
    print("init_gmix successful. Features:")
    for(i in 1:n_feats){
      print(Parm$Features$Names[i])
    }
  }
  return(list(Parm=Parm))
}

#
#
#
#
#