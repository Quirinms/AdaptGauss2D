\name{computeDotDensity}
\alias{computeDotDensity}
\title{computeDotDensity}
\description{
computeDotDensity
}
\usage{
computeDotDensity(Data, Means, Covariances, Weights)
}
\arguments{
\item{Data}{[1:n, 1:2] Numeric matrix with n observations and 2 feature
dimensions.}
\item{Means}{List of l two-dimensional numeric vectors carrying mean vector.}
\item{Covariances}{List of l 2x2 numeric matrices carrying covariance matrices.}
\item{Weights}{[1:l] Numeric vector with l weights for l Gaussian Mixture Model
components.}
}
\value{
\item{DotDensity}{[1:n] Numeric vector with density values for each observation
from Data.}
}
\author{Quirin Stier}
\keyword{GMM}
\concept{gaussian mixture model}
\keyword{EM}
\concept{expectation maximization}
\keyword{multimodal}
\keyword{mixture}
\keyword{Baggenstoss}
