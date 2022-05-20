\name{covariance2AxesAngle}
\alias{covariance2AxesAngle}
\title{covariance2AxesAngle}
\description{
Transforms two numbers (variances of independent directions) and an angle to an
ellipsoid matrix.
}
\usage{
covariance2AxesAngle(Covariances)
}
\arguments{
\item{Covariances}{[1:2] Variance from two independent directions (data axes or
more simple: standard deviation of data from first and second dimensions
respectively).}
}
\value{
\item{MainAxes}{[1:3] Numeric vector}
}
\author{Quirin Stier}
\keyword{GMM}
\concept{gaussian mixture model}
\keyword{EM}
\concept{expectation maximization}
\keyword{multimodal}
\keyword{mixture}
\keyword{Baggenstoss}
