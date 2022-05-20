EstimateNumberOfModes=function(Data,MaxModeNo=7,SampleSize=4000,...){
  if(!requireNamespace("FCPS")){
    warning("EstimateNumberOfModes: FCPS package is missing")
    return(2)
  }
  if(!requireNamespace("parallelDist")){
    warning("EstimateNumberOfModes: parallelDist package is missing")
    return(2)
  }
  #sollte war dann nicht gehen
  # if(!requireNamespace("fastcluster")){
  #   #for ward algorithm
  #   warning("EstimateNumberOfModes: fastcluster package is missing")
  #   return(2)
  # }
  if(!is.matrix(Data)){
    warning("EstimateNumberOfModes: Data is expected to be a matrix. Calling as.matrix()...")
    Data=as.matrix(Data)
  }
  if(mode(Data)!="numeric"){
    warning("EstimateNumberOfModes: Data is expected to be numeric. transfering to numeric...")
    mode(Data)=="numeric"
  }
  if(length(SampleSize)!=1){
    warning("EstimateNumberOfModes: SampleSize is expected to be numeric of length one, trying to take first element..")
    SampleSize=SampleSize[1]
  }

  n=nrow(Data)
  if(n>SampleSize){
    Data=Data[sample(1:n,SampleSize,replace = F),]
  }
  DD=parallelDist::parallelDist(Data)
  hc <- hclust(DD, method = "ward.D2")
  clsm <- matrix(data = 0, nrow = dim(Data)[1],ncol = MaxModeNo)
  for (i in 2:(MaxModeNo+1)) {
    clsm[,i-1] <- cutree(hc,i)
  }
  colnames(clsm) = 1:dim(clsm)[2]
  ModesNoVec=c()
  for(i in 1:20){
    res=FCPS::ClusterNoEstimation(Data = Data, ClsMatrix = clsm,MaxClusterNo =  MaxModeNo,ClusterIndex = "ssi",...)
    ModesNoVec[i]=as.vector(res$ClusterNo)
  }
  ModesNo=round(mean(ModesNoVec,na.rm=TRUE),0)

  return(list(ModesNo=ModesNo,ListObj=res))
}
