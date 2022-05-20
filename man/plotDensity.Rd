\name{plotDensity}
\alias{plotDensity}
\title{Plot density of two-dimensional data}
\description{
Plot density of two-dimensional data.
}
\usage{
plotDensity(Data, Cls, CurrGauss, Colors,

GridDensity,Means, Covariances,Weights,

MainAxesAngle,XKernel, YKernel, Shapes, ShapeText,

AxNames = c("X1", "X2"),ShowAxis = FALSE, ShowEllipsoids = TRUE,

ShowScatter = FALSE,ShowGaussNr = FALSE, Source = "D")
}
\arguments{
\item{Data}{[1:n, 1:2] Numeric matrix with n observations and 2 features.}
\item{Cls}{[1:n] Numerical vector of size n containing the classes of each observation.}
\item{CurrGauss}{Integer number indicating the currently highlighted gaussian
component.}
\item{Colors}{...}
\item{GridDensity}{[1:x, 1:x] Numeric matrix with PDF for a GMM with one or more
Gaussians on a Grid defined by the vectors XKernel and YKernel.}
\item{Means}{List with l [1:2] numerical vector defining the means of the l GMM
components.}
\item{Covariances}{List with l [1:2, 1:2] numerical matrices defining the covariance
matrices of the l GMM components.}
\item{Weights}{[1:l] Numerical vector with weights for each GMM component.}
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
\item{ShowAxis}{Boolean determining ...}
\item{ShowEllipsoids}{Boolean determining ....}
\item{ShowScatter}{Boolean determining ....}
\item{ShowGaussNr}{Boolean determining if the shown Gaussian components are
enumerated (TRUE) or not (FALSE).}
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
