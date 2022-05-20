covariance2AxesAngle = function(Covariances){
  #
  # DESCRIPTION
  # Computes the PCA information of a list of 2x2 covariance matrices. These
  # computations yield the PCA vectors, weights and the angle of the larger
  # principal component axis with respect to the cartesian coordinate system.
  #
  # INPUT
  # Cov[1:2, 1:2]    Numerical matrix containing the covariance matrix.
  #
  # OUTPUT
  # PCAinfo    List
  #
  # DETAILS
  # The length of the principial component axes of a covariance matrix are
  # computed with the singular value decomposition. The angle is calculated
  # here with a cosine formula and the dot product. The angle here is measured
  # as the angle between the first (larger) principal component and the
  # canoncial vector (1,0). Starting from the cartesian coordinate (1,0)
  # the computed angle degrees range from 0 to 360 in counter clock direction.
  #
  # Author QS 2021

  if(missing(Covariances)){
    message("Parameter Cov is missing. Returning.")
    return()
  }else{
    if(!is.matrix(Covariances)){
      message("Parameter Cov is not of type matrix. Returning.")
      return()
    }
  }
  # The SVD (singular value decomposition) computes the main axes of a matrix
  # Here: we want the main axes of a covariance matrix (should always be positive definit and symmetric)
  # Then we can always obtain an ellipsoid
  # Note: We always obtain the largest ax as first element (order of values is decreasing)
  MySVD = svd(Covariances)                                   # Compute singular value decomposition for Princ. Component Axes
  PCWeight1 = MySVD$d[1]                                  # Extract 1st singular value (first PCA component)
  PCWeight2 = MySVD$d[2]                                  # Extract 2nd singular value
  PCVector1 = MySVD$u[,1]                                 # Extract 1st PCA component vector
  PCVector2 = MySVD$u[,2]                                 # Extract 2nd PCA component vector
  NormPC1 = norm(PCVector1, type = "2")
  TopCircle1    = acos(sum(PCVector1 * c(0,1))/NormPC1)*(180/pi) # See if 1st PCA is on the upper part of the cartesian coord. sys.
  BottomCircle1 = acos(sum(PCVector1 * c(0,-1))/NormPC1)*(180/pi)
  Angle1        = acos(sum(PCVector1 * c(1,0))/NormPC1)*(180/pi)
  if(BottomCircle1<TopCircle1){
    Angle1 = 360 - Angle1                                   # This would be the angle for the lower part
  }
  NormPC2 = norm(PCVector2, type = "2")
  TopCircle2    = acos(sum(PCVector2 * c(0,1))/NormPC2)*(180/pi) # See if 1st PCA is on the upper part of the cartesian coord. sys.
  BottomCircle2 = acos(sum(PCVector2 * c(0,-1))/NormPC2)*(180/pi)
  Angle2        = acos(sum(PCVector2 * c(1,0))/NormPC2)*(180/pi)
  if(BottomCircle2<TopCircle2){
    Angle2 = 360 - Angle2                                   # This would be the angle for the lower part
  }
  return(c(sqrt(PCWeight1), sqrt(PCWeight2), Angle1, Angle2))
}
#
#
#
#
#
