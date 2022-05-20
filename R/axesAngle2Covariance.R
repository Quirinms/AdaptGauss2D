axesAngle2Covariance = function(Variance2D, Angle){
  #
  # DESCRIPTION
  # Computes the covariance matrix given the desired standard deviation
  # of the main axes of the ellipsoid an the angle between the main ax
  # of the ellipsoid and the standard cartesian system. The angle must be
  # measured from the canonical unit vector (1,0).
  #
  # INPUT
  # Axes[1:2]          Variances from two independent directions (data axes or
  #                    more simple: standard deviations of data from first and
  #                    second dimensions respectively).
  # Angle              Numerical value defining the angle between the main ax of
  #                    an ellipsoid to the standard cartesian system.
  #
  # OUTPUT
  # CovarianceMatrix2D[1:2, 1:2]    Numerical matrix containing the covariance matrix.
  #
  # Author QS 2021
  if(missing(Variance2D)){
    message("Parameter Variance2D is missing. Returning.")
    return()
  }else{
    if(!is.vector(Variance2D)){
      message("Parameter Variance2D is not of type vector. Returning.")
      return()
    }else if(length(Variance2D) != 2){
      message("Parameter Variance2D must be a vector of dimension 2. Returning.")
      return()
    }
  }
  if(missing(Angle)){
    message("Parameter Angle is missing. Returning.")
    return()
  }else if(!is.numeric(Angle)){
    message("Parameter Angle must be a numeric value. Returning.")
    return()
  }

  # Draw ellipse via two axes X,Y with length x,y
  # Rotate ellipse by angle Angle
  # This results in Covariance matrix
  #InternalVariance   = c(max(Variance2D), min(Variance2D))
  Trafo              = 180/pi                                         # Degree to Radiant
  RotMatEntries      = c(cos(Angle/Trafo), -sin(Angle/Trafo),         # Rotation matrix
                         sin(Angle/Trafo),  cos(Angle/Trafo))
  RotationMatrix     = matrix(RotMatEntries, ncol = 2, byrow = T)
  CovarianceMatrix2D = RotationMatrix %*% (diag(Variance2D)**2) %*% t(RotationMatrix) # Axes length is sqrt of eigenvalues

  #varX1 = majorAxis² * cos(phi)² + minorAxis² * sin(phi)²
  #varX2 = majorAxis² * sin(phi)² + minorAxis² * cos(phi)²
  #cov12 = (majorAxis² - minorAxis²) * sin(phi) * cos(phi)

  #Trafo = pi/180
  #CovarianceMatrix2D = matrix(0,2,2)
  #CovarianceMatrix2D[1,1] = Variance2D[1]**2 * cos(Angle*Trafo)**2 + Variance2D[2]**2 * sin(Angle*Trafo)**2
  #CovarianceMatrix2D[1,2] = (Variance2D[1]**2 - Variance2D[2]**2) * sin(Angle*Trafo) * cos(Angle*Trafo)
  #CovarianceMatrix2D[2,1] = CovarianceMatrix2D[1,2]
  #CovarianceMatrix2D[2,2] = Variance2D[1]**2 * sin(Angle*Trafo)**2 + Variance2D[2]**2 * cos(Angle*Trafo)**2
  return(round(CovarianceMatrix2D, digits = 4))
}
#
#
#
#
#
