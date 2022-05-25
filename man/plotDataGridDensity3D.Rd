\name{plotDataGridDensity3D}
\alias{plotDataGridDensity3D}
\title{Plot Density of two-dimensional Data in three dimensions}
\description{
Visualize the density of a two-dimensional dataset in three dimensions.
}
\usage{
plotDataGridDensity3D(Data, XKernel, YKernel, ContinuousDataPDE, EmpiricDataPDE,
Colors, Cls, Camera = NULL, ShowScatter = TRUE, AxNames = c("X", "Y"),
Source = "F1", Debug = F)
}
\arguments{
\item{Data}{[1:n, 1:2] Numeric matrix with n observations and 2 features.}
\item{XKernel}{[1:x] Numeric vector defining domain of x axis.}
\item{YKernel}{[1:x] Numeric vector defining domain of y axis.}
\item{ContinuousDataPDE}{[1:x, 1:x] Numeric matrix with density estimation of
Data computed as grid on the domain defined by XKernel and YKernel in order to
plot a continuous density plot.}
\item{EmpiricDataPDE}{[1:n] Numeric vector with density estimation of Data
defined for each datapoint within the Data.}
\item{Colors}{[1:l] Character vector with l color names or more.}
\item{Cls}{[1:n] Numeric vector with class assignment for n observations.}
\item{Camera}{List of attributes concerning the camera angle for plotly visualizations.}
\item{ShowScatter}{Boolean (Default=TRUE). T: Show Data Scatter with empirical
estimated density. F: Do not show nothing.}
\item{AxNames}{[1:2] Character vector with ax names for plot description.}
\item{Source}{Character indicating plot source. Important attribute for
plotly in shiny in order to keep control of specific panels.}
\item{Debug}{Boolean (Default=FALSE). TRUE: Show developer information and warnings
in terminal. FALSE: Show nothing.}
}
\value{
\item{plotOut}{Plotly object containing plot for direct visualization.}
}
\author{Quirin Stier}
\examples{

\dontrun{
data(EngyTime)
Data = EngyTime
densityMap = MASS::kde2d(Data[,1], Data[,2], n = 100)
z = sapply(1:nrow(Data), function(i){
  densityMap$z[
    which.min(abs(densityMap$x - Data[i,1])),
    which.min(abs(densityMap$y - Data[i,2]))
    ]
})
dataDensity = z
Cls = rep(1, nrow(Data))
DataSample = sample(1:nrow(Data), min(nrow(Data), 5000))
Colors = c("red", "blue")
Plot_3D_DataDensity(Data, DataSample, Colors, Cls, densityMap, dataDensity,
Show3DPoints = T, Source = "F1")
}
}
\keyword{GMM}
\concept{gaussian mixture model}
\keyword{EM}
\concept{expectation maximization}
\keyword{multimodal}
\keyword{mixture}
