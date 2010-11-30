\name{mse}
\Rdversion{1.1}
\alias{mse}
\alias{mse.default}
\alias{mse.matrix}
\alias{mse.data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean Squared Error
}
\description{
Mean squared error between \code{sim} and \code{obs}, in the squared units of \code{sim} and \code{obs}, with treatment of missing values.
}
\usage{
mse(sim, obs, ...)

\method{mse}{default}(sim, obs, na.rm=TRUE, ...)

\method{mse}{data.frame}(sim, obs, na.rm=TRUE, ...)

\method{mse}{matrix}(sim, obs, na.rm=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sim}{
numeric, zoo, matrix or data.frame with simulated values
}
  \item{obs}{
numeric, zoo, matrix or data.frame with observed values
}
  \item{na.rm}{
a logical value indicating whether 'NA' should be stripped before the computation proceeds. \cr
When an 'NA' value is found at the i-th position in \code{obs} \bold{OR} \code{sim}, the i-th value of \code{obs} \bold{AND} \code{sim} are removed before the computation.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Mean squared error between \code{sim} and \code{obs}. \cr

If \code{sim} and \code{obs} are matrixes, the returned value is a vector, with the mean squared error between each column of \code{sim} and \code{obs}.
}
\references{
\cite{Yapo P. O., Gupta H. V., Sorooshian S., 1996. Automatic calibration of conceptual rainfall-runoff models: sensitivity to calibration data. Journal of Hydrology. v181 i1-4. 23-48}
}
\author{
Mauricio Zambrano Bigiarini <mauricio.zambrano@ing.unitn.it>
}
\note{
\code{obs} and \code{sim} has to have the same length/dimension \cr

The missing values in \code{obs} and \code{sim} are removed before the computation proceeds, and only those positions with non-missing values in \code{obs} and \code{sim} are considered in the computation \cr
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mae}}, \code{\link{me}}
}
\examples{
obs <- 1:10
sim <- 1:10
mse(sim, obs)

obs <- 1:10
sim <- 2:11
mse(sim, obs)

##################
# Loading daily streamflows of the Ega River (Spain), from 1961 to 1970
require(zoo)
data(EgaEnEstellaQts)
obs <- EgaEnEstellaQts

# Generating a simulated daily time series, initially equal to the observed series
sim <- obs 

# Computing the mean squared error for the "best" case
mse(sim=sim, obs=obs)

# Randomly changing the first 2000 elements of 'sim', by using a normal distribution 
# with mean 10 and standard deviation equal to 1 (default of 'rnorm').
sim[1:2000] <- obs[1:2000] + rnorm(2000, mean=10)

# Computing the new mean squared error
mse(sim=sim, obs=obs)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }