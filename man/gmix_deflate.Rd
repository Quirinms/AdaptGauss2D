\name{gmix_deflate}
\alias{gmix_deflate}
\title{EM Baggenstoss: gmix_deflate}
\description{
Subroutine to update Gaussian mixture (5 Operations):
2. Pruning modes (gmix_deflate.R)
The method gmix_deflate.R is killing weak modes (a mode is another name for one
of the L mixture components). A weak mode is found by testing the correct mode
entry of weight vector wts to see if it falls below a threshold. Only one mode
per run can be killed.
See the matlab documentation for more information.
}
\usage{
gmix_deflate(Parm,MinWeight1,MinWeightAll,Verbose=0)
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
\item{MinWeight1}{Numerical value: First mode below this threshold is killed.}
\item{MinWeightAll}{Numerical value: All modes below this threshold are
killed.}
\item{Verbose}{Optional: Don't print comments. Verbose=1 print them. Default=0.}
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
