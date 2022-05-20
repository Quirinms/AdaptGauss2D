\name{axesAngle2Covariance}
\alias{axesAngle2Covariance}
\title{Axes and Angle to Ellipsoid Matrix}
\description{
Transforms two numbers (variances of independent directions) and an angle to an
ellipsoid matrix.
}
\usage{
axesAngle2Covariance(Variance2D, Angle)
}
\arguments{
\item{Variance2D}{[1:2] Variance from two independent directions (data axes or
more simple: standard deviation of data from first and second dimensions
respectively).}
\item{Angle}{Numerical value defining the angle between the main ax of an
ellipsoid to the standard cartesian system.}
}
\value{
\item{CovarianceMatrix2D}{[1:2, 1:2] Numeric matrix containing the covariance
matrix.}
}
\author{Quirin Stier}
\keyword{GMM}
\concept{gaussian mixture model}
\keyword{EM}
\concept{expectation maximization}
\keyword{multimodal}
\keyword{mixture}
\keyword{Baggenstoss}
