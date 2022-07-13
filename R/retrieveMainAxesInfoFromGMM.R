retrieveMainAxesInfoFromGMM = function(Covariance, MainAxesAngle){
  #
  # DESCRIPTION
  #
  #
  # INPUT
  # Covariances                    List with l [1:2, 1:2] numerical matrices
  #                                defining the covariance matrices of the l GMM
  #                                components.
  # MainAxesAngle[1:4]             Numeric vector with 1st and 2nd main axes
  #                                of a 2D ellipsoid and the respective angles
  #                                measured to the first unit vector c(0,1).  #
  # OUTPUT
  # PC1A    Coordinate of the first main ax in the first dimension of two.
  # PC1B    Coordinate of the first main ax in the second dimension of two.
  # PC2A    Coordinate of the second main ax in the first dimension of two.
  # PC2B    Coordinate of the second main ax in the second dimension of two.
  #
  # Author QS 2021
  MySVD = svd(Covariance)                                   # Compute singular value decomposition for Princ. Component Axes
  SD1 = MySVD$d[1]*MySVD$u[,1]                                 # Extract 1st PCA component vector
  SD2 = MySVD$d[2]*MySVD$u[,2]                                 # Extract 2nd PCA component vector
  NormSD1   = norm(SD1, type = "2")
  TopCircle1    = acos(sum(SD1 * c(0,1))/NormSD1)*(180/pi)  # See if 1st PCA is on the upper part of the cartesian coord. sys.
  BottomCircle1 = acos(sum(SD1 * c(0,-1))/NormSD1)*(180/pi) # See if 1st PCA is on the lower part of the cartesian coord. sys.
  Angle1 = acos(sum(SD1 * c(1,0))/NormSD1)*(180/pi)
  if(BottomCircle1<TopCircle1){
    Angle1 = 360 - Angle1                                   # This would be the angle for the lower part
  }
  if(round(abs(Angle1-MainAxesAngle[3])) > 5 & MainAxesAngle[3] != 360){
    SD1 = -SD1
  }
  NormSD2   = norm(SD2, type = "2")
  TopCircle2    = acos(sum(SD2 * c(0,1))/NormSD2)*(180/pi) # See if 1st PCA is on the upper part of the cartesian coord. sys.
  BottomCircle2 = acos(sum(SD2 * c(0,-1))/NormSD2)*(180/pi)
  Angle2 = acos(sum(SD2 * c(1,0))/NormSD2)*(180/pi)
  if(BottomCircle2<TopCircle2){
    Angle2 = 360 - Angle2                                   # This would be the angle for the lower part
  }
  if(abs(Angle2-((MainAxesAngle[3]+90)%%360)) > 5){
    SD2 = -SD2
  }
  PC1A = SD1[1]; PC1B = SD1[2]; PC2A = SD2[1]; PC2B = SD2[2] # Eigenvector components
  return(list("PC1A" = PC1A,
              "PC1B" = PC1B,
              "PC2A" = PC2A,
              "PC2B" = PC2B))
}
#
#
#
#
#
