\name{computeRMS}
\alias{computeRMS}
\title{computeRMS}
\description{
computeRMS
}
\usage{
computeRMS(Data, Means, Covariances, Weights, EmpiricDataPDE)
}
\arguments{
\item{Data}{[1:n, 1:2] Numeric matrix with n observations and 2 feature
dimensions.}
\item{Means}{List of l two-dimensional numeric vectors carrying mean vector.}
\item{Covariances}{List of l 2x2 numeric matrices carrying covariance matrices.}
\item{Weights}{[1:l] Numeric vector with l weights for l Gaussian Mixture Model
components.}
\item{EmpiricDataPDE}{[1:n] Numeric vector with density estimation of Data
defined for each datapoint within the Data.}
}
\value{
\item{Mu}{List of l two-dimensional numeric vectors carrying mean vector.}
\item{Sigma}{List of l 2x2 numeric matrices carrying covariance matrices.}
\item{Lambda}{[1:l] Numeric vector with l weights for l Gaussian Mixture Model
components.}
\item{MainAxesAngle}{List of length l with numeric vectors [1:3] with main axes and the angle
between the first dimension of a 2D cartesian coordinate system the first main
ax for all Gaussian Mixture Model components.}
}
\author{Quirin Stier}
\keyword{GMM}
\concept{gaussian mixture model}
\keyword{EM}
\concept{expectation maximization}
\keyword{multimodal}
\keyword{mixture}
\keyword{Baggenstoss}
