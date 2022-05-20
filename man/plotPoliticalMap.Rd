\name{plotPoliticalMap}
\alias{plotPoliticalMap}
\title{Political Map of Gaussian Mixture Model}
\description{
Draws a political map of the classification of a Gaussian Mixture Model.
}
\usage{
plotPoliticalMap(Data, CurrGauss,Means, Covariances, Weights,

MainAxesAngle,Colors, Cls,Shapes, ShapeText, AxNames = c("X1", "X2"),

ShowAxis = FALSE, ShowEllipsoids = TRUE,

 ShowScatter = FALSE, Source = "D")
}
\arguments{
\item{Data}{[1:n, 1:2] Numeric matrix with n observations and 2 features.}
\item{CurrGauss}{...}
\item{Means}{List with l [1:2] numerical vector defining the means of the l GMM
components.}
\item{Covariances}{List with l [1:2, 1:2] numerical matrices defining the
covariance matrices of the l GMM components.}
\item{Weights}{[1:l] Numerical vector with weights for each GMM component.}
\item{MainAxesAngle}{List of numeric vectors with 1st and 2nd main axes of a 2D
ellipsoid and the respective angles measured to the first unit vector c(0,1).}
\item{Colors}{[1:l] Character vector with l color names.}
\item{Cls}{[1:n] Numeric vector with classification for n observations.}
\item{Shapes}{List of List with 4 attributes (type, fillcolor, opacity, path)
for a shape for plotting. Here it is used for plotting an ellipsoid.}
\item{ShapeText}{[1:l, 1:3] Numeric matrix with l means and two entries for the
two-dimensional coordinates and one entry for the number of the Gaussian
component.}
\item{AxNames}{[1:2] Character vector with ax names for plot description.}
\item{ShowAxis}{Boolean indicating if data scatter is shown or not.}
\item{ShowEllipsoids}{Boolean indicating if data scatter is shown or not.}
\item{ShowScatter}{Boolean indicating if data scatter is shown or not.}

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
\keyword{Baggenstoss}
