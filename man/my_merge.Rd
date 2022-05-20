\name{my_merge}
\alias{my_merge}
\title{EM Baggenstoss: my_merge}
\description{
Merges two gmix modes.
}
\usage{
my_merge(Weight1, Mean1, Covariance1, Weight2, Mean2, Covariance2)
}
\arguments{
\item{Weight1}{Numerical value with weight of first mode.}
\item{Mean1}{Numerical vector with mean of first mode (d-dimensional where d
denotes the features dimension).}
\item{Covariance1}{Numerical matrix with covariance matrix of first mode (dxd matrix
where d denotes the feature's dimension).}
\item{Weight2}{Numerical value with weight of second mode.}
\item{Mean2}{Numerical vector with mean of second mode  (d-dimensional where
d denotes the features dimension).}
\item{Covariance2}{Numerical matrix with covariance matrix of second mode (dxd
matrix where d denotes the feature's dimension).}
}
\value{
List with
\item{Weight}{Numerical value with weight of merged mode.}
\item{Mean}{Numerical vector with mean of merged mode (d-dimensional where
d denotes the features dimension).}
\item{Covariance}{Numerical matrix with covariance matrix of merged mode (dxd
matrix where d denotes the feature's dimension).}
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
