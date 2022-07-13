\name{EstimateNumberOfModes}
\alias{EstimateNumberOfModes}
\title{
Estimates the number of modes
}
\description{
Estimates the number of modes for two-dimensional gaussian distributed data using [Thrun/Stier, 2021].
}
\usage{
EstimateNumberOfModes(Data, MaxModeNo = 7, SampleSize = 4000, ...)
}
\arguments{
  \item{Data}{
[1:n,1:d] numerical matrix of a dataset to be clustered. It consists of n cases of d-dimensional data points. Every case has d attributes, variables or features.
}
  \item{MaxModeNo}{
    Highest number of modes to be checked
}
  \item{SampleSize}{
size of sample to be taken that is tested for the number of modes
}
  \item{\dots}{
Optional further arguments to ClusterNoEstimation
}
}
\details{
An empirical approach is taken internally by using a combination of ward and ssi index.
}
\value{
LIST of
\item{ModesNo }{scalar, numerical value describing the number of clusters}
\item{ListObj }{resulting object of function  ClusterNoEstimation of package FCPS}
}
\references{
[Thrun/Stier, 2021]  Thrun, M. C., & Stier, Q.: Fundamental Clustering Algorithms Suite SoftwareX, Vol. 13(C), pp. 100642, doi 10.1016/j.softx.2020.100642, 2021.
}
\author{
Michael Thrun
}
\seealso{
\link[FCPS]{ClusterNoEstimation}
}
\examples{
data(EngyTime)
EstimateNumberOfModes(EngyTime$Data)$ModesNo
}
\keyword{modes}
\concept{number of modes}
