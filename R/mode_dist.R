mode_dist = function(Mean1, Mean2, Covariance1, Covariance2){
  # V = mode_dist(Mean1, Mean2, Covariance1, Covariance2)
  # V$parm
  # 
  # DESCRIPTION
  # Distance (actually closeness) measure between 2 modes
  # with means Mean1, Mean2 and cholesky covariance matrices Covariance1,
  # Covariance2 returns negative closeness measure (for identical modes, d = 0)
  # 
  # 
  # INPUT
  # Mean1[1:d]              Numerical vector with d feature dimensions
  # Mean2[1:d]              Numerical vector with d feature dimensions
  # Covariance1[1:d,1:d]    Numerical matrix dxd (d features)
  # Covariance2[1:d,1:d]    Numerical matrix dxd (d features)
  # 
  # 
  # OUTPUT
  # Dist     Dist
  # 
  # 
  #  Dr. Paul M. Baggenstoss
  #  Naval Undersea Warfare Center
  #  Newport, RI
  #  p.m.baggenstoss@ieee.org
  # 
  # in /dbt/EMforGauss/
  # Author QS 2021
  
  #Mean1 = as.vector(t(m1))
  #m2 = as.vector(t(m2))
  #v1 = reshape(v1,dim,dim);
  #v2 = reshape(v2,dim,dim);
  
  dim     = length(Mean1)
  npoints = 2*(2*dim+1)
  x       = matrix(0, nrow = npoints, ncol = dim)
  point   = 1
  
  for(imode in 0:1){
    # select the right means and variances
    if(imode == 0){
      tr = Covariance1
      mean = Mean1
    }else{
      tr = Covariance2
      mean = Mean2
    }
    # Get the eigenvectors (right singular vectors (columns of V) 
    res_svd = svd(tr)
    S = res_svd$d
    V = res_svd$v
    # Add the center point 
    x[point,] = mean
    point = point + 1
    # loop over eigenvectors 
    for(ev in 1:dim){
      x[point,] = mean + S[ev] * V[,ev]
      point = point + 1
      x[point,] = mean - S[ev] * V[,ev]
      point = point + 1
    }
  }
  V    = lqr_eval(x,Mean1,Covariance1)
  out1 = V$LogPDF
  V    = lqr_eval(x,Mean2,Covariance2)
  out2 = V$LogPDF
  
  a = sum(out1[((npoints/2)+1):npoints])    # p1(x2)
  b = sum(out2[1:(npoints/2)])              # p2(x1)
  c = sum(out2[((npoints/2)+1):npoints])    # p2(x2)
  d = sum(out1[1:(npoints/2)])              # p1(x1)
  Dist = a+b-c-d
  return(list(Dist=Dist))
}
#
#
#
#
#
