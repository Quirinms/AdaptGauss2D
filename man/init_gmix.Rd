\name{init_gmix}
\alias{init_gmix}
\title{EM Baggenstoss: init_gmix}
\description{
Initialization of GMM parameters for training on data with gmix_trainscript.R
}
\usage{
init_gmix(Data, Nmodes, MinStd, Names=NULL, RandInit=1, Verbose=0)
}
\arguments{
\item{Data}{[1:n,1:d] Numerical matrix with n samples and d feature dimensions.}
\item{Nmodes}{[1:L] Number of modes to initialize.}
\item{MinStd}{Numerical vector with covariance constraints for each feature.
Preventing the covariance matrix to become singular.}
\item{Names}{[1:d] String vector with feature names.}
\item{RandInit}{Optional: Use randomized initialization. Random=1,
Nonrandom=0.}
\item{Verbose}{Optional: Print Output. 1=Yes, 0=No. Default=1.}
}
\value{
List with one named element parm, which is a nested list
\item{Parm}{Nested list with parameters for GMM.
Features carrying permanent values.
Features$name [1:d] String vector with feature names.
Features$min_std [1:NMODE] Vector of covariance constraints.
Modes carrying modifyable values.
Modes$cholesky_covar [d*NMODE, d] Numerical matrix with NMODE many square
matrices stacked vertically with the covariance matrix.
Modes$mean [1:NMODE, d] Numerical matrix with nmode different means and d
feature dimensions.
Modes$weight [1, 1:NMODE] Numerical matrix with weights for each mean.}
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
