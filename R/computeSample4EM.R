computeSample4EM = function(N, BigDataMaximum = 0.6){
  # DESCRIPTION
  #
  # INPUT
  # N                 Integer meaning the number of observations.
  # BigDataMaximum    Numeric between 0 and 1 meaning the percent of data to use.
  #
  # OUTPUT
  # dataSampleIdx[1:l]    Numeric vector with sample of data for visualization.
  #
  # Author: QMS 13.04.2022
  # Use only a part of the data for EM
  # For Big Data (>= 10k), use 60% of the data for EM
  #N = nrow(Data)
  tmpNumSample1   = round(max(min(N, 2000), BigDataMaximum*N), digits = 0)
  dataSampleIdx   = sample(1:N, tmpNumSample1)
  return(dataSampleIdx)
}
