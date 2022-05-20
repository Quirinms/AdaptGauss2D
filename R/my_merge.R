my_merge = function(Weight1, Mean1, Covariance1, Weight2, Mean2, Covariance2){
  # V = my_merge(Weight1, Mean1, Covariance1, Weight2, Mean2, Covariance2)
  # V$Weight
  # V$Mean
  # V$Covariance
  # 
  # DESCRIPTION
  # Merges two gmix modes
  # 
  # INPUT
  # Weight1                Numerical value with weight of first mode.
  # Mean1[1:d]             Numerical vector with mean of first mode
  #                        (d-dimensional where d denotes the features
  #                        dimension).
  # Covariance1[1:d,1d]    Numerical matrix with covariance matrix of first mode
  #                        (dxd matrix where d denotes the feature's dimension).
  # Weight2                Numerical value with weight of second mode.
  # Mean2[1:d]             Numerical vector with mean of second mode
  #                        (d-dimensional where d denotes the features
  #                        dimension).
  # Covariance2[1:d,1d]    Numerical matrix with covariance matrix of second
  #                        mode (dxd matrix where d denotes the feature's
  #                        dimension).
  # 
  # OUTPUT
  # Weight                Numerical value with weight of merged mode.
  # Mean[1:d]             Numerical vector with mean of merged mode 
  #                       (d-dimensional where d denotes the features
  #                       dimension).
  # Covariance[1:d,1d]    Numerical matrix with covariance matrix of merged mode
  #                       (dxd matrix where d denotes the feature's dimension).
  # 
  # 
  #  Dr. Paul M. Baggenstoss
  #  Naval Undersea Warfare Center
  #  Newport, RI
  #  p.m.baggenstoss@ieee.org
  # 
  # in /dbt/EMforGauss/
  # Author QS 2021
  DIM     = length(Mean1)
  Weight  = Weight1 + Weight2
  nor     = (Weight1 + Weight2)
  Weight1 = Weight1/nor
  Weight2 = Weight2/nor
  # The first step is to determine the central mean by weighted average
  Mean = Weight1 * Mean1 + Weight2 * Mean2
  # the ROWS of v times the square root of DIM (v is the the QR of covariance)
  # can be considered data samples about the means of each distribution.
  Covariance1 = Covariance1 * sqrt(DIM)
  Covariance2 = Covariance2 * sqrt(DIM)
  # They need to be re-referenced to the new center
  for(i in 1:DIM){
    Covariance1[i,] = Covariance1[i,] + (Mean1 - Mean)
    Covariance2[i,] = Covariance2[i,] + (Mean2 - Mean)
  }
  matTemp    = rbind(sqrt(Weight1)*Covariance1, sqrt(Weight2)*Covariance2)   # form a weighted augmented matrix
  res_qr     = qr(matTemp)
  tmpvar     = qr.R(res_qr)
  Covariance = tmpvar * (1/sqrt(DIM))
  return(list(Weight=Weight, Mean=Mean, Covariance=Covariance))
}
#
#
#
#
#