\name{mode_dist}
\alias{mode_dist}
\title{EM Baggenstoss: mode_dist}
\description{
Distance (actually closeness) measure between 2 modes with means m1, m2 and
cholesky covariance matrices v1, v2. Returns negative closeness measure
(for identical modes, d = 0).
}
\usage{
mode_dist(Mean1, Mean2, Covariance1, Covariance2)
}
\arguments{
\item{Mean1}{Numerical vector with mean of first mode (d-dimensional where d
denotes the features dimension).}
\item{Mean2}{Numerical vector with mean of second mode  (d-dimensional where
d denotes the features dimension).}
\item{Covariance1}{Numerical matrix with covariance matrix of first mode (dxd matrix
where d denotes the feature's dimension).}
\item{Covariance2}{Numerical matrix with covariance matrix of second mode (dxd
matrix where d denotes the feature's dimension).}
}
\value{
List with one element:
\item{Dist}{Numerical value with distance of two modes.}
}
\author{Quirin Stier}
\references{
Baggenstoss, Paul M., and T. E. Luginbuhl.: An EM algorithm for joint model estimation. IEEE International Conference on Acoustics, Speech, and Signal Processing. Proceedings. ICASSP99 (Cat. No.99CH36258), Phoenix, AZ, USA, 1999, pp. 1825-1828 vol.4, IEEE, doi:10.1109/ICASSP.1999.758276, 1999.
}
\keyword{GMM}
\keyword{gaussian mixture model}
\keyword{EM}
\keyword{expectation maximization}
\keyword{multimodal}
\keyword{mixture}
\keyword{Baggenstoss}
