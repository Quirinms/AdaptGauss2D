\name{lqr_eval}
\alias{lqr_eval}
\title{EM Baggenstoss: lqr_eval}
\description{
Calculates the Gaussian logPDF given Mean and Covariances.
}
\usage{
lqr_eval(Data, ModeMean, CholeskyCovar)
}
\arguments{
\item{Data}{[1:N,1:d] Numerical matrix, N = number of cases, d = dimension.}
\item{ModeMean}{Numerical vector with one d-dim mode in GMM.}
\item{CholeskyCovar}{Numerical matrix with covariance matrix of d-dim mode
(dxd matrix where d denotes the feature's dimension).}
}
\value{
List with
\item{logPDF}{Numerical vector with log Gaussian PDF.}
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
