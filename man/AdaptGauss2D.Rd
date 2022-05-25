\name{AdaptDunes}
\alias{AdaptDunes}
\title{Adapt a two-dimensional Gaussian Mixture Model}
\description{
GUI for manually fitting a multivariate (2 dimensional) gaussian mixture to a
dataset.
}
\usage{
AdaptDunes(Data, Means = NULL, CovarianceMatrices = NULL, Weights = NULL,
Cls = NULL, Debug = F, dbt = F, WorkingDirectory = getwd())
}
\arguments{
\item{Data}{[1:n, 1:2] Numeric matrix with n observations and 2 feature
dimensions.}
\item{Means}{List of l two-dimensional numeric vectors carrying mean vector.}
\item{CovarianceMatrices}{List of l 2x2 numeric matrices carrying covariance
matrices.}
\item{Weights}{[1:l] Numeric vector with l weights for l Gaussian Mixture Model
components.}
\item{Cls}{[1:n] Integer vector with n entries for classes assigned to each
data point in Data.}
\item{Debug}{Boolean (Default=FALSE). T: Show developer information and warnings
in terminal. F: Show nothing.}
\item{dbt}{Boolean (Default=FALSE). T: Provide more complex functionality to the
dbt. F: Hide unnecessary details to keep it simple for a broad audience
especially from the non computer science community.}
\item{WorkingDirectory}{Character String (Default=getwd()). Indicates the
working directory to which to save the setting obtained with AdaptDunes to.}
}
\value{
List with multiple elements of different data structures:
\item{Means}{List of l numerical 2D vectors. Each vector representing the
expectation values for one Gaussian component.}
\item{CovarianceMatrices}{List of l numerical 2x2 matrices.}
\item{Weights}{[1:l] Numerical vector of size l containing the weight for each
Gaussian component within the mixture.}
\item{PrincipalComponentAxis}{List of l numerical 2D vectors. Each matrix
representing the prinicpal component axis of the covariance matrix of a Gaussian
component.}
\item{Angle}{List of l integers. Each integer represents the angle of the first
principal component ax of the covariance matrix.}
\item{Cls}{[1:n] Numerical vector of size n containing the classes of each observation.}
}
\details{
The AdaptDunes app is a GUI made for an intuitive approach for two-dimensional
gaussian mixture modeling (GMM). There is also an approach for the one-dimensional
case, called AdaptGauss.
The output of the app can be used as input. Be aware, that the principal
component axis and the angle are not required, since this information is
within the covariance matrices and is only provided in the output as a quick
access of the results.
The essence of AdaptDunes lies in the coarse manual fitting as a first step and
a fine automatized fitting with an EM algorithm as the second step.
That way the user can communicate a rough expectation of the model which is
then optimized. In other words, the user provides a rough desired initialization
and the automaton of AdaptDunes finds a local optimum for that.
AdaptDunes yields an initialized solution for a GMM as proposal and start for
the user.
}
\author{Quirin Stier}
\examples{
\dontrun{
data(EngyTime)
# Call shiny application
Data=EngyTime$Data

Res = AdaptDunes(Data, Means = NULL, CovarianceMatrices = NULL, Weights = NULL,
Cls = NULL, Debug = F)
}
}
\keyword{GMM}
\concept{gaussian mixture model}
\keyword{EM}
\concept{expectation maximization}
\keyword{multimodal}
\keyword{mixture}
\keyword{Baggenstoss}
