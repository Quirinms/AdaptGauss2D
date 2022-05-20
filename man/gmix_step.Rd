\name{gmix_step}
\alias{gmix_step}
\title{EM Baggenstoss: gmix_step}
\description{
Subroutine to update Gaussian mixture (5 Operations):
1. E-M algorithm (gmix_step)
E-M algorithm (expectation maximization algorithm)
See the matlab documentation for more information.
}
\usage{
gmix_step(Parm, Data, Bias=0, DataWTS=NULL)
}
\arguments{
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
\item{Data}{[1:n,1:d] Numerical matrix with normalized data. N samples with DIM 
feature dimensions.}
\item{Bias}{Optional: Default=0. Binary value: Covariance constraint method.
Choose: 1=BIAS, 0=CONSTRAINT}
\item{DataWTS}{Optional: Numerical Vector, which allows individually weighting
input data. Default=matrix(1,N,1).}
}
\value{
List with two elements:
\item{Parm}{Nested list with parameters for GMM.
                    Parm$features carrying permanent values.
                    
                    Parm$features$name [1:d]        String vector with
                                                    feature names.
                    Parm$features$min_std [1:NMODE] Vector of covariance
                                                    constraints.
                                                    
                    Parm$modes carrying modifyable values.
                    Parm$modes$cholesky_covar [d*NMODE, d]     Numerical
                         matrix with NMODE many square matrices stacked 
                         vertically with the covariance matrix.
                    Parm$modes$mean [1:NMODE, d]    Numerical matrix with nmode
                         different means and d feature dimensions.
                    Parm$modes$weight [1, 1:NMODE]    Numerical matrix with
                         weights for each mean.}
\item{Q}{Numerical value: Total log-likelihood output (weighted by data_wts).}
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
