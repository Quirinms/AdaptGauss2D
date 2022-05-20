\name{ad_gmix_trainscript}
\alias{ad_gmix_trainscript}
\title{Interface to Baggenstoss'es EM-Algorithm to calculate Gaussian Mixture Model}
\description{
Interface to Baggenstoss'es EM-Algorithm to calculate Gaussian Mixture Model for
given data.
}
\usage{
ad_gmix_trainscript(Parm, Data, Nit, SamplesPerMode=NULL, Bias=0, Maxclose=NULL,
Addmodes=1, KurtosisThreshold=1.0, Verbose=0)
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
\item{Nit}{Numerical value: max number of iterations.}
\item{SamplesPerMode}{Optional: Numerical value:
Samples-per-mode (minimum for pruning). Default: 4*d.}
\item{Bias}{Optional: Binary value: Covariance constraint method.
Choose: 1=BIAS, 0=CONSTRAINT. Default=0.}
\item{Maxclose}{Optional: Numerical value: Maximum mode closeness. Use
more negative values to promote mode consolidation. Use higher values for larger
dimension (Default: Suggestion: -2 * d).}
\item{Addmodes}{Optional: 0 or 1 value: If set to 1, will use kurt.m to split
modes. Default = 1.}
\item{KurtosisThreshold}{Optional: Default=1.0. Numerical value: Kurtosis and Skew
threshold for mode splitting. Should be about 1.0.  Higher values (i.e. 1.2)
will make mode splitting less likely.}
\item{Verbose}{Optional: Optional: 0, 1 or 2. 0 = no output, 1 = messages, 2 =
messages + plot. Default ==0. Print some outputs.}
}
\value{
List with multiple elements of different data structures:
\item{EMmean}{[1:l,1:d] Numerical matrix carrying L means row wise.}
\item{EMCov}{[1:(l*d),1:d] Numerical matrix with all covariance matrices appended
beneath each other.}
\item{EMalpha}{[1:l] Vector with wights for all L modes.}
\item{EMParams}{Nested list with parameters for GMM. EMParams is the adapted or
trained input variable Parm, see input argument Parm.}
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
