\name{ad_interfaceEM}
\alias{ad_interfaceEM}
\title{ad_interfaceEM}
\description{
Transforms two numbers (variances of independent directions) and an angle to an
ellipsoid matrix.
}
\usage{
ad_interfaceEM(Data, Mu, Sigma, Lambda, Nmodes = 1, Addmodes = 1, Nit = 1,
Verbose = 0)
}
\arguments{
\item{Data}{[1:n, 1:2] Numeric matrix with n observations and 2 feature
dimensions.}
\item{Mu}{List of l two-dimensional numeric vectors carrying mean vector.}
\item{Sigma}{List of l 2x2 numeric matrices carrying covariance matrices.}
\item{Lambda}{[1:l] Numeric vector with l weights for l Gaussian Mixture Model
components.}
\item{Nmodes}{Number of modes to initialize if no (Mu, Sigma, Lambda) is given.}
\item{Nit}{Numerical value: max number of iterations.}
\item{Addmodes}{Optional: Addmodes = 1: Allow EM to add modes. Addmodes = 0: do
not allow EM to add modes.}
\item{Verbose}{Optional: 0, 1 or 2. 0 = no output, 1 = messages, 2 =
messages + plot. Default ==0. Print some outputs.}
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
