decimal2CompassDegree = function(DecimalDegree){
  #
  # DESCRIPTION
  #
  #
  # INPUT
  #
  #
  # OUTPUT
  #
  #
  # DETAILS
  #
  # Author QS 2021
  if((DecimalDegree >= 0 ) & (DecimalDegree <= 90)){
    CompassDegree = 90 - DecimalDegree
  }else if((DecimalDegree >= 90 ) & (DecimalDegree <= 180)){
    CompassDegree = 270 + 180 - DecimalDegree
  }else if((DecimalDegree >= 180) & (DecimalDegree <= 270)){
    CompassDegree = 180 + 90 - (DecimalDegree - 180)
  }else{
    CompassDegree = 180 - (DecimalDegree - 270)
  }
  return(CompassDegree)
}
#
#
#
#
#
