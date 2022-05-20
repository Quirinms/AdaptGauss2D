computeSample4Visualization = function(N){
  # DESCRIPTION
  #
  # INPUT
  # N                 Integer meaning the number of observations.
  #
  # OUTPUT
  # dataSampleIdx[1:l]    Numeric vector with sample of data for visualization.
  #
  # Author: QMS 13.04.2022

  # Use only a part of the data for visualization purposes
  # Use only 2000 datapoints
  tmpNumSample  = min(N, 2000)#round(max(min(nrow(Data), 2000), 0.1*nrow(Data)), digits = 0)
  dataSampleIdx = sample(1:N, tmpNumSample)
  return(dataSampleIdx)
}
