\name{br2}
\Rdversion{1.1}
\alias{br2}
\alias{br2.default}
\alias{br2.matrix}
\alias{br2.data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
br2
}
\description{
Coefficient of determination (r2) multiplied by the slope of the regression line between \code{sim} and \code{obs}, with treatment of missing values.
}
\usage{
br2(sim, obs, ...)

\method{br2}{default}(sim, obs, na.rm=TRUE, ...)

\method{br2}{data.frame}(sim, obs, na.rm=TRUE, ...)

\method{br2}{matrix}(sim, obs, na.rm=TRUE, ...)
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
\details{
\deqn{br2 = \mid b \mid R2, \mid b \mid \leq 1; br2= \mid b \mid R2, \mid b \mid > 1 }{%
br2 = abs(b)*R2, b <= 1; br2= abs(b)/R2, abs(b) > 1 } 

A model that systematically over or under-predicts all the time will still result in "good" \code{r2} (close to 1), even if all predctions were wrong (Krause et al., 2005).
The \code{br2} coefficient allows accounting for the discrepancy in the magnitude of two signals (depicted by 'b') as well as their dynamics (depicted by 'r2')
}
\value{
br2 between \code{sim} and \code{obs}. \cr

If \code{sim} and \code{obs} are matrixes, the returned value is a vector, with the br2 between each column of \code{sim} and \code{obs}.
}
\references{
\cite{Krause, P., Boyle, D. P., and Base, F.: Comparison of different efficiency criteria for hydrological model assessment, Adv. Geosci., 5, 89-97, 2005}

}
\author{
Mauricio Zambrano Bigiarini <mauricio.zambrano@ing.unitn.it>
}
\note{
\code{obs} and \code{sim} has to have the same length/dimension \cr

The missing values in \code{obs} and \code{sim} are removed before the computation proceeds, and only those positions with non-missing values in \code{obs} and \code{sim} are considered in the computation \cr

The slope 'b' is computed as the coefficient of the linear regression between \code{sim} and \code{obs}, forcing the intercept be equal to zero.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gof}}, \code{\link{cor}}, \code{\link{lm}}
}
\examples{
obs <- 1:10
sim <- 1:10
br2(sim, obs)

obs <- 1:10
sim <- 2:11
br2(sim, obs)

##################
# Loading daily streamflows of the Ega River (Spain), from 1961 to 1970
require(zoo)
data(EgaEnEstellaQts)
obs <- EgaEnEstellaQts

# Generating a simulated daily time series, initially equal to the observed series
sim <- obs 

# Computing 'br2' for the "best" (unattainable) case
br2(sim=sim, obs=obs)

# Randomly changing the first 2000 elements of 'sim', by using a normal distribution 
# with mean 10 and standard deviation equal to 1 (default of 'rnorm').
sim[1:2000] <- obs[1:2000] + rnorm(2000, mean=10)

# Computing the new  'br2'
br2(sim=sim, obs=obs)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }