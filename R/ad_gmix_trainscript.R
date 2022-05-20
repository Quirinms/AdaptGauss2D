ad_gmix_trainscript = function(Parm, Data, Nit, SamplesPerMode=NULL, Bias=0,
                               Maxclose=NULL, Addmodes=1, KurtosisThreshold=1.0,
                               Verbose=0){
  # V = gmix_trainscript(Parm, Data, Nit, SamplesPerMode=NULL, Bias=0,
  #                      Maxclose=NULL, Addmodes=1, KurtosisThreshold=1.0, Verbose=0)
  # V$EMmean
  # V$EMCov
  # V$EMalpha
  # V$EMParams
  #
  # DESCRIPTION
  # Interface to Baggenstoss'es EM-Algorithm to calculate Gaussian Mixture Model
  # for given Data
  #
  #
  # INPUT
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
  # Data[1:n,1:d]    Numerical matrix with normalized data. n samples with d
  #                  feature dimensions.
  # Nit              Numerical value: max number of iterations.
  #
  # OPTIONAL
  # SamplesPerMode       Numerical value:Samples-per-mode (minimum for pruning).
  #                      Default: 4*d
  # Bias                 Binary value: Covariance constraint method.
  #                      Choose: 1=BIAS, 0=CONSTRAINT.
  #                      Default = 0.
  # Maxclose             Optional: Numerical value: Maximum mode closeness. Use
  #                      more negative values to promote mode consolidation. Use
  #                      higher values for larger dimension.
  #                      Default: -2 * DIM.
  #                      Suggestion: -0.5 times DIM.
  # Addmodes             0 or 1 value: If set to 1, will use kurt.m to split modes.
  #                      Default = 1.
  # KurtosisThreshold    Numerical value: Kurtosis and Skew threshold for mode
  #                      splitting. Should be about 1.0.  Higher values
  #                      (i.e. 1.2) will make mode splitting less likely.
  # Verbose              Optional: 0, 1 or 2.
  #                      0 = no output, 1 = messages, 2 = messages + plot.
  #                      Default ==0. Print some outputs.
  #
  #
  #
  # OUTPUT
  # EMmean[1:l,1:d]       Numerical matrix with l row vectors containing the d
  #                       dimensional means.
  # EMCov[1:d,1:d,1:l]    Numerical matrix  with l many square matrices stacked
  #                       vertically with the covariance matrix.
  # EMalpha[1:l]          Numerical vector with weights for each mean.
  # EMParams              Nested list with parameters for GMM.
  #                       Parm$Features carrying permanent values.
  #                       Parm$Features$Names [1:d]      String vector with
  #                                                      feature names.
  #                       Parm$Features$MinStd [1:l]     Vector of covariance
  #                                                       constraints.
  #                       Parm$Modes list of modifyable GMM parameters for all modes.
  #                       Parm$Modes$Covariance [1:d*l, 1:d]     Numerical matrix
  #                            with l many square matrices stacked vertically with
  #                            the covariance matrix.
  #                       Parm$Modes$Mean [1:l, 1:d]                 Numerical matrix
  #                           with l row vectors containing the d dimensional means.
  #                       Parm$modes$Weight [1:l]                    Numerical vector
  #                            Modes weights for each mean.
  #
  #
  # Dr. Paul M. Baggenstoss
  # Naval Undersea Warfare Center
  # Newport, RI
  # p.m.baggenstoss@ieee.org
  #
  # in /dbt/EMforGauss/
  # Author QS 2021
  #----------------------------------------------------------------------------#
  # Check if required variables are there with correct data structure
  if(is.null(Parm)){
    stop("Input parameter Parm is not given. Returning.")
  }
  if(is.null(Data)){
    stop("Input parameter Data is not given. Returning.")
  }
  if(is.null(Nit)){
    stop("Input parameter Nit is not given. Returning.")
  }
  if(!is.matrix(Data)){
    stop("Data is not a matrix!")
  }
  if(!is.numeric(Nit)){
    stop("Nit must be a numerical value!")
  }
  if(!is.list(Parm)){
    stop("Parm must be a list!")
  }
  #----------------------------------------------------------------------------#
  # Check if parameter Parm has all required entries
  if(is.null(Parm$Features)){
    stop("Parm must contain sublist Parm$Features!")
  }
  if(is.null(Parm$Features$Names)){
    stop("Parm must contain String vector Parm$Features$Names!")
  }
  if(is.null(Parm$Features$MinStd)){
    stop("Parm must contain numerical vector Parm$Features$MinStd!")
  }
  if(is.null(Parm$Modes$Covariance)){
    stop("Parm must contain numerical matrix Parm$Modes$Covariance!")
  }
  if(is.null(Parm$Modes$Mean)){
    stop("Parm must contain numerical matrix Parm$Modes$Mean!")
  }
  if(is.null(Parm$Modes$Weight)){
    stop("Parm must contain numerical vector Parm$Modes$Weight!")
  }
  #----------------------------------------------------------------------------#
  # Check if all entries of Parm have correct data structure
  if(!is.list(Parm$Features)){
    stop("Parm$Features must be a list!")
  }
  if(!is.character(Parm$Features$Names)){
    stop("Parm$Features$Names must be a String vector (is.character)!")
  }
  if(!is.double(Parm$Features$MinStd)){
    stop("Parm$Features$MinStd must be a numerical vector!")
  }
  if(!is.double(Parm$Modes$Covariance)){
    stop("Parm$Modes$Covariance must be a numerical matrix!")
  }
  if(!is.double(Parm$Modes$Mean)){
    stop("Parm$Modes$Mean must be a numerical matrix!")
  }
  if(!is.double(Parm$Modes$Weight)){
    stop("Parm$Modes$Covariance must be a numerical vector!")
  }
  #----------------------------------------------------------------------------#
  # Check if dimensionalities check out
  if(dim(Parm$Modes$Covariance)[2] != dim(Data)[2]){
    stop("Feature dimensions of covariance matrix and Data differ but must be
         equal!")
  }
  # Number of modes must check out, both in number of weights and number of
  # means
  if(length(Parm$Modes$Weight) != dim(Parm$Modes$Mean)[1]){
    stop("Number of weights of mode (see Parm$Modes$Weight) and number of means
         dim(Parm$Modes$Mean)[1] differ but must be equal!")
  }
  #----------------------------------------------------------------------------#
  wts   = Parm$Modes$Weight    #
  nmode = length(wts)          # determine the number of modes
  N     = dim(Data)[1]
  DIM   = dim(Data)[2]
  if(is.null(SamplesPerMode)){
    SamplesPerMode = 4*DIM
  }
  if(is.null(Maxclose)){
    Maxclose = -2 * DIM
  }
  stop_crit    = .0001*N;  # How little an improvement is needed to continue
  train_period = 6;        # How often to check for close modes and try mode splitting
  if(nmode < 3){
    imerge = 0
    isplit = round(train_period/2)
  }else{
    imerge = round(train_period/2)
    isplit = 0
  }
  nq = 0                      # reset convergence counter
  if((DIM*N*Nit)<(1E6)){    # Small problem do it without tacho
    for(iter in 1:Nit){       # Nit max number of iterations
      V    = gmix_step(Parm, Data, Bias) # Q Total log-likelihood output
      Parm = V$Parm
      Q    = V$Q
      #if((Verbose == 2)&(iter > 5)){
      #  if(Addmodes != 0){ # adding / gausses is allowed
      #    #subplot
      #    #hold on; plot(iter,-log(-Q),'.'); title('costs'); hold off; axis tight;
      #    #grid on;
      #    #subplot(1,2,2);
      #    #hold on; plot(iter,nmode,'.'); title(['Nr. of Gaussians: ',num2str(nmode)]); hold off; axis tight;
      #    #grid on;
      #    #drawnow;
      #  }else{
      #    #hold on; plot(iter,abs(-log(-Q)),'.'); title('EM costs'); hold off; axis tight;
      #    #grid on;
      #  }
      #}else{
      #  if(Verbose == 1){
      #    print(paste("Q(", iter, ")=", Q, "Nr. of Gaussians=", nmode, sep=""))
      #  }
      #}
      if(iter == 1){
        best = Q
      }
      wts   = Parm$Modes$Weight
      nmode = length(wts)
      if(Addmodes){ # evenually reduce number of gaussians ALU 2019
        V    = gmix_deflate(Parm = Parm,
                            MinWeight1 = min(SamplesPerMode/N,.5),
                            MinWeightAll = 1.5/N,
                            Verbose = Verbose)
        Parm = V$Parm
      }
      if(iter < (Nit-train_period)){
        if((iter > (2*train_period)) & (iter%%train_period)==imerge){
          V    = gmix_merge(Parm = Parm,
                            MaxCloseness = Maxclose,
                            Verbose = Verbose)
          Parm = V$Parm
        }
        if(Addmodes & ((iter%%train_period) == isplit)){
          V    = gmix_kurt(Parm = Parm,
                           Data = Data,
                           Threshold = KurtosisThreshold,
                           Debug = 0)
          Parm = V$Parm
        }
      }
      wts   = Parm$Modes$Weight
      if(length(wts) != nmode){      # If there has been a change in modes, allow more time to converge
        nq = 0
      }
      nmode = length(wts)
      if((Q - best) > stop_crit){    # has there been an improvement?
        nq = 0                       # reset counter if there's an improvement
      }else{
        # Increment "no improvement" counter. If no improvement
        # has been seen in train_period steps, terminate
        nq = nq+1
        if(nq > train_period){
          break
        }
      }
      if(Q > best){
        best = Q
      }
    }
  }else{ # (DIM*N*Nit)> 1E6  large problem use tacho
    UpdateTacho = 10                     # when to update tacho
    Pname ='EM '                         # namen fuer die tacho Ausgabe
    #ShowTacho = Verbose > 0

    withProgress(message = 'Computing EM', value = 0, {
      for(iter in 1:Nit){
        V    = gmix_step(Parm, Data, Bias)
        Parm = V$Parm
        Q    = V$Q
        LogQ = abs(-log(-Q))

        incProgress(1/Nit, detail = paste("Doing step", iter))

        #if ShowTacho &(mod(iter,UpdateTacho)==0) ;
        #waitbar(iter/Nit,Tacho,[Pname,': ',num2str(iter),' /', num2str(Nit),' Q=',num2str(round(LogQ))]);
        #end; % update tacho
        if((Verbose ==2) & ( iter>5)){
          if(Addmodes != 0){ # adding / gausses is allowed
            #plot(iter, -log(-Q))
            #subplot(1,2,1);
            #hold on; plot(iter,LogQ,'.'); title('costs'); hold off; axis tight;
            #grid on;
            #subplot(1,2,2);
            #hold on; plot(iter,nmode,'.'); title(['Nr. of Gaussians: ',num2str(nmode)]); hold off; axis tight;
            #grid on;
            #drawnow;
          }else{
            plot(iter, -log(-Q))
            #hold on; plot(iter,LogQ,'.'); title('EM costs'); hold off; axis tight;
            #grid on;
          }
        }else{
          if(Verbose ==1){
            print(paste("Q(", iter, ")=", Q, "Nr. of Gaussians=", nmode, sep=""))
            #print("")
          }
        }
        if(iter==1){
          best = Q
        }
        wts   = Parm$Modes$Weight
        nmode = length(wts)

        if(Addmodes){ # evenually reduce number of gaussians ALU 2019
          V    = gmix_deflate(Parm, min(SamplesPerMode/N,.5), 1.5/N,Verbose)
          Parm = V$Parm
        }
        if(iter < Nit-train_period){
          if((iter > 2*train_period) & ((iter%%train_period)==imerge)){
            V    = gmix_merge(Parm, Maxclose,Verbose)
            Parm = V$Parm
          }
          if(Addmodes & ((iter%%train_period) == isplit)){
            V    = gmix_kurt(Parm, Data, KurtosisThreshold, 0)
            Parm = V$Parm
          }
        }
        wts   = Parm$Modes$Weight
        # If there has been a change in modes, allow more time to converge
        if(length(wts) != nmode){
          nq=0
        }
        nmode=length(wts)
        # has there been an improvement?
        if( Q - best > stop_crit ){
          # reset counter if there's an improvement
          nq = 0
        }else{
          #Increment "no improvement" counter. If no improvement
          # has been seen in train_period steps, terminate
          nq = nq+1
          if(nq > train_period){
            if(Verbose==TRUE){
              print("Break condition reached. Finishing training early.")
            }
            break
          }
        }
        if(Q > best){
          best=Q
        }
        #if(ShowTacho){
        #  close(pb)
        #}
      }
    })
  }
  colnames(Parm$Modes$Mean) = NULL
  return(list("EMmean"   = Parm$Modes$Mean,
              "EMCov"    = Parm$Modes$Covariance,
              "EMalpha"  = Parm$Modes$Weight,
              "EMParams" = Parm))
}
#
#
#
#
#
