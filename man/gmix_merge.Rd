\name{gmix_merge}
\alias{gmix_merge}
\title{EM Baggenstoss: gmix_merge}
\description{
Subroutine to update Gaussian mixture (5 Operations):
3. Merging modes (gmix_merge.R)
This method creates a single mode from two nearly identical modes.
The closeness of two modes is determined by mode_dist.R.
See the matlab documentation for more information.
}
\usage{
gmix_merge(Parm,MaxCloseness, Verbose=0)
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
\item{MaxCloseness}{Optional: Numerical value: Maximum mode closeness. Use
more negative values to promote mode consolidation. Use higher values for larger
dimension (suggest -0.5 times DIM). Default=-2 * DIM.}
\item{Verbose}{Optional: Print some outputs. Default=0.}
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
