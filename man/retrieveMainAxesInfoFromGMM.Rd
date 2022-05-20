\name{retrieveMainAxesInfoFromGMM}
\alias{retrieveMainAxesInfoFromGMM}
\title{retrieveMainAxesInfoFromGMM}
\description{
retrieveMainAxesInfoFromGMM
}
\usage{
retrieveMainAxesInfoFromGMM(Covariance, MainAxesAngle)
}
\arguments{
\item{Covariance}{List with l [1:2, 1:2] numerical matrices defining the covariance
matrices of the l GMM components.}
\item{MainAxesAngle}{List of numeric vectors with 1st and 2nd main axes of a 2D
ellipsoid and the respective angles measured to the first unit vector c(0,1).}
}
\value{
\item{PC1A}{Coordinate of the first main ax in the first dimension of two.}
\item{PC1B}{Coordinate of the first main ax in the second dimension of two.}
\item{PC2A}{Coordinate of the second main ax in the first dimension of two.}
\item{PC2B}{Coordinate of the second main ax in the second dimension of two.}
}
\author{Quirin Stier}
\keyword{GMM}
\concept{gaussian mixture model}
\keyword{EM}
\concept{expectation maximization}
\keyword{multimodal}
\keyword{mixture}
\keyword{Baggenstoss}
