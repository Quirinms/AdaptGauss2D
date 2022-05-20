\name{gmix_kurt}
\alias{gmix_kurt}
\title{EM Baggenstoss: gmix_kurt}
\description{
Subroutine to update Gaussian mixture (5 Operations):
4. Splitting modes (gmix_kurt.R)
A split of a mode can be determined with the weighted kurtosis for each mode for
one dimensional problems.
For multidimensional problems the skew is used in addition to the kurtosis.
See the matlab documentation for the explicit mathematical formulas and remarks.
}
\usage{
gmix_kurt(Parm,Data,Threshold,Debug=0)
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
\item{Threshold}{Optional: Numerical value: . Default=1.0.}
\item{Debug}{Optional: Print some outputs. Default=0.}
}
\value{
List with one named element parm, which is a nested list
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
