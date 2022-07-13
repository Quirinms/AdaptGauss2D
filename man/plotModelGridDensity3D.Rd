\name{plotModelGridDensity3D}
\alias{plotModelGridDensity3D}
\title{Plot Gaussian Mixture Model density of two-dimensional data in three
dimensions}
\description{
Compares the continuous density plot of a dataset with the density distribution
defined by one or more gaussians (GMM; Mean, Cov, Weights).
The datapoints from the dataset have the density assigned by the GMM component
with the the highest probability (class assignment must be given and is not
computed here).
}
\usage{
plotModelGridDensity3D(Data, XKernel, YKernel, ContinuousDataPDE, GridDensity,
Means, Covariances, Weights, MainAxesAngle,
Colors, Cls,
Camera = NULL, AxNames = c("X", "Y"), Source = "F1", Debug = FALSE)
}
\arguments{
\item{Data}{[1:n, 1:2] Numeric matrix with n observations and 2 features.}
\item{XKernel}{[1:x] Numeric vector defining domain of x axis.}
\item{YKernel}{[1:x] Numeric vector defining domain of y axis.}
\item{ContinuousDataPDE}{[1:x, 1:x] Numeric matrix with density estimation of
Data computed as grid on the domain defined by XKernel and YKernel in order to
plot a continuous density plot.}
\item{GridDensity}{[1:x, 1:x] Numeric matrix with PDF for a GMM with one or more
Gaussians on a Grid defined by the vectors XKernel and YKernel.}
\item{Means}{List with l [1:2] numerical vector defining the means of the l GMM
components.}
\item{Covariances}{List with l [1:2, 1:2] numerical matrices defining the covariance
matrices of the l GMM components.}
\item{Weights}{[1:l] Numerical vector with weights for each GMM component.}
\item{MainAxesAngle}{List of numeric vectors with 1st and 2nd main axes of a 2D
ellipsoid and the respective angles measured to the first unit vector c(0,1).}
\item{Colors}{[1:l] Character vector with l color names or more.}
\item{Cls}{[1:n] Numerical vector of size n containing the classes of each observation.}
\item{Camera}{List of attributes concerning the camera angle for plotly visualizations.}
\item{AxNames}{[1:2] Character vector with ax names for plot description.}
\item{Source}{(Default = "F1"). Character indicating plot source. Important attribute for
plotly in shiny in order to keep control of specific panels.}
\item{Debug}{Boolean (Default=FALSE). TRUE: Show developer information and warnings
in terminal. FALSE: Show nothing.}
}
\value{
\item{plotOut}{Plotly object containing plot for direct visualization.}
}
\details{
Uses a density estimation (discretized computation on grid) to display
a continuous density plot on domain defined by XKernel and YKernel of the
whole empircal data distribution.
The density of the GMM components is computed via an multivariate gaussian
formula.
The datapoints take the density of the GMM components with highest
probability to take the specific datapoint. This computation must be done
outside of this routine (give the Cls/Class/Labels as parameter).
Parameter Camera is important for shiny to let the view fixed during
continuous updates of the interactive visualization (Shiny App). The Source
is important to separate between multiple visualizations in the Shiny App.
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
