\name{lqr_evp}
\alias{lqr_evp}
\title{EM Baggenstoss: lqr_evp}
\description{
Helper function for creating classed based on Bayes decision with Expectation
Maximization.
Returns either total log-likelihood of Gaussian mixtures or separate log
likelihoods for all modes.
}
\usage{
lqr_evp(Parm, Data, Flag)
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
\item{Flag}{Boolean flag. flag==1, returns log mode PDF in columns, that is the
output of lqr_eval only, not including mixing weights for each column.
Flag=0: Computes the total log-PDF for normalized data in a single output column
(i.e. the log of the weighted sum of the rows of the flag=1 output) plus the
Jacobian (it outputs the  PDF value with respect to the raw unnormalized data by
taking into account the jacobian of the normalization operation). In this case
outputs ModePDFs and PDF are derived.}
}
\value{
List with
\item{lg}{[1:n] Numerical vector with LogPDF}
\item{ModePDFs}{[1:n,1:L] Numerical matrix with PDFs for all single modes.}
\item{PDF}{[1:n] Numerical vector with combined Gaussian as weighthed sum of
PDF.}
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
