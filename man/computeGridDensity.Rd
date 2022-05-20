\name{computeGridDensity}
\alias{computeGridDensity}
\title{computeGridDensity}
\description{
computeGridDensity
}
\usage{
computeGridDensity(XKernel, YKernel, Means, Covariances, Weights)
}
\arguments{
\item{XKernel}{[1:x] Numeric vector defining domain of x axis.}
\item{YKernel}{[1:x] Numeric vector defining domain of y axis.}
\item{Means}{List of l two-dimensional numeric vectors carrying mean vector.}
\item{Covariances}{List of l 2x2 numeric matrices carrying covariance matrices.}
\item{Weights}{[1:l] Numeric vector with l weights for l Gaussian Mixture Model
components.}
}
\value{
\item{GridDensity}{[1:x, 1:x] Numeric matrix with PDF for a GMM with one or more
Gaussians on a Grid defined by the vectors XKernel and YKernel.}
}
\author{Quirin Stier}
\keyword{GMM}
\concept{gaussian mixture model}
\keyword{EM}
\concept{expectation maximization}
\keyword{multimodal}
\keyword{mixture}
\keyword{Baggenstoss}
