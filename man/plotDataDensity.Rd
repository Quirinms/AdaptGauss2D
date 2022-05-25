\name{plotDataDensity}
\alias{plotDataDensity}
\title{Plot density of two-dimensional data}
\description{
Plot density of two-dimensional data.
}
\usage{
plotDataDensity(Data, Cls, CurrGauss, Colors,
ContinuousDataPDE,Means, Covariances, MainAxesAngle,
XKernel, YKernel, Shapes, ShapeText, AxNames = c("X1", "X2"),
ShowAxis = FALSE, ShowEllipsoids = TRUE, ShowScatter = FALSE,
ShowGaussNr = FALSE, Source = "D")
}
\arguments{
\item{Data}{[1:n, 1:2] Numeric matrix with n observations and 2 features.}
\item{Cls}{[1:n] Numerical vector of size n containing the classes of each observation.}
\item{CurrGauss}{Integer indicating the current selected Gaussian
component.}
\item{Colors}{[1:l] Character vector with l color names or more.}
\item{ContinuousDataPDE}{[1:x, 1:x] Numeric matrix with density estimation of
Data computed as grid on the domain defined by XKernel and YKernel in order to
plot a continuous density plot.}
\item{Means}{List with l [1:2] numerical vector defining the means of the l GMM
components.}
\item{Covariances}{List with l [1:2, 1:2] numerical matrices defining the covariance
matrices of the l GMM components.}
\item{MainAxesAngle}{List of numeric vectors with 1st and 2nd main axes of a 2D
ellipsoid and the respective angles measured to the first unit vector c(0,1).}
\item{XKernel}{[1:x] Numeric vector defining domain of x axis.}
\item{YKernel}{[1:x] Numeric vector defining domain of y axis.}
\item{Shapes}{List of List with 4 attributes (type, fillcolor, opacity, path)
for a shape for plotting. Here it is used for plotting an ellipsoid.}
\item{ShapeText}{[1:l, 1:3] Numeric matrix with l means and two entries for the
two-dimensional coordinates and one entry for the number of the Gaussian
component.}
\item{AxNames}{[1:2] Character vector with ax names for plot description.}
\item{ShowAxis}{Boolean indicating to show models axis or not.}
\item{ShowEllipsoids}{Boolean indicating to show models ellipsoids or not.}
\item{ShowScatter}{Boolean indicating if data scatter is shown or not.}
\item{ShowGaussNr}{Boolean indicating to enumerate model components.}
\item{Source}{Character indicating plot source (Default = "D"). Important attribute for
plotly in shiny in order to keep control of specific panels.}
}
\value{
\item{plotOut}{Plotly object containing plot for direct visualization.}
}
\author{Quirin Stier}
\keyword{GMM}
\concept{gaussian mixture model}
\keyword{EM}
\concept{expectation maximization}
\keyword{multimodal}
\keyword{mixture}
