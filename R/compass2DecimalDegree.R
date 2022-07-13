compass2DecimalDegree = function(CompassDegree){
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
  if((CompassDegree >= 0 ) & (CompassDegree <= 90)){
    DecimalDegree = 90 - CompassDegree
  }else if((CompassDegree >= 90 ) & (CompassDegree <= 180)){
    DecimalDegree = 270 + 180 - CompassDegree
  }else if((CompassDegree >= 180) & (CompassDegree <= 270)){
    DecimalDegree = 180 + 90 - (CompassDegree - 180)
  }else{
    DecimalDegree = 180 - (CompassDegree - 270)
  }
  return(DecimalDegree)
}
#
#
#
#
#
