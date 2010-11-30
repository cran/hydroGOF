\name{gof}
\Rdversion{1.1}
\alias{gof}
\alias{gof.default}
\alias{gof.matrix}
\alias{gof.data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Numerical Goodness of Fit
}
\description{
Numerical goodness of fit between \code{sim} and \code{obs}, with treatment of missing values. Several performance indexes for comparing two vectors, matrixes or data.frames
}
\usage{
gof(sim, obs, ...)

\method{gof}{default}(sim, obs, na.rm=TRUE, do.spearman=FALSE, do.pbfdc=FALSE, digits=2,...)
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
  \item{do.spearman}{
logical. Indicates if the Spearman correlation has to be computed. The default is FALSE.
}
  \item{do.pbfdc}{
logical. Indicates if the Percent Bias in the Slope of the Midsegment of the Flow Duration Curve (\code{\link{pbiasfdc}}) has to be computed. The default is FALSE.
}
  \item{digits}{
decimal places used for rounding the goodness-of-fit indexes.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
%%\details{
%%}
\value{
\item{me}{Mean Error}
\item{mae}{Mean Absolute Error}
\item{mse}{Mean Squared Error}
\item{rmse}{Root Mean Square Error}
\item{nrmse}{Normalized Root Mean Square Error  ( -100\% <= nrms <= 100\% )}
\item{PBIAS}{Percent Bias }
\item{pbfdc}{PBIAS in the slope of the midsegment of the Flow Duration Curve}
\item{RSR}{Ratio of RMSE to the Standard Deviation of the Observations, RSR = rms / sd(obs). ( 0 <= RSR <= +Inf )}
\item{rSD}{Ratio of Standard Deviations, rSD = sd(sim) / sd(obs)}
\item{NSeff}{Nash-Sutcliffe Efficiency ( -Inf <= NSeff <= 1 )}
\item{NSeff1}{Modified Nash-Sutcliffe Efficiency}
\item{rNSeff}{Relative Nash-Sutcliffe Efficiency }
\item{d}{Index of Agreement ( 0 <= d <= 1 )}
\item{d1}{Modified Index of Agreement}
\item{rd}{Relative Index of Agreement}
\item{cp}{Persistence Index ( 0 <= PI <= 1 ) }
\item{r}{Pearson Correlation coefficient ( -1 <= r <= 1 )}
\item{r.Spearman}{Spearman Correlation coefficient ( -1 <= r.Spearman <= 1 ) }
\item{R2}{ Coefficient of Determination ( 0 <= R2 <= 1 ). \cr
          Gives the proportion of the variance of one variable that is predictable from the other variable}
\item{bR2}{R2 multiplied by the coefficient of the regression line between \code{sim} and \code{obs} \cr ( 0 <= bR2 <= 1 )}
%% ...
}
\references{
\cite{Legates, D. R., and G. J. McCabe Jr. (1999), Evaluating the Use of "Goodness-of-Fit" Measures in Hydrologic and Hydroclimatic Model Validation, Water Resour. Res., 35(1), 233--241} \cr

\cite{Krause P., Boyle D.P., and B\"{a}se F., Comparison of different efficiency criteria for hydrological model assessment, Advances in Geosciences 5 (2005), pp. 89--97} \cr

\cite{Moriasi, D.N., Arnold, J.G., Van Liew, M.W., Bingner, R.L., Harmel, R.D., Veith, T.L. 2007. Model evaluation guidelines for systematic quantification of accuracy in watershed simulations \cr Transactions of the ASABE. 50(3):885-900} \cr

\cite{Boyle, D. P., H. V. Gupta, and S. Sorooshian (2000), Toward Improved Calibration of Hydrologic Models: Combining the Strengths of Manual and Automatic Methods, Water Resour. Res., 36(12), 3663--3674} \cr

\cite{Kitanidis, P. K., and R. L. Bras (1980), Real-Time Forecasting With a Conceptual Hydrologic Model 2. Applications and Results, Water Resour. Res., 16(6), 1034--1044} \cr

\cite{J.E. Nash and J.V. Sutcliffe, River flow forecasting through conceptual models. Part 1: a discussion of principles, J. Hydrol. 10 (1970), pp. 282--290} \cr

\cite{Yapo P. O., Gupta H. V., Sorooshian S., 1996. Automatic calibration of conceptual rainfall-runoff models: sensitivity to calibration data. Journal of Hydrology. v181 i1-4. 23--48} \cr

\cite{Yilmaz, K. K., H. V. Gupta, and T. Wagener  (2008), A process-based diagnostic approach to model evaluation: Application to the NWS distributed hydrologic model, Water Resour. Res., 44, W09417, doi:10.1029/2007WR006716}
}
\author{
Mauricio Zambrano Bigiarini <mauricio.zambrano@ing.unitn.it>
}
\note{
\code{obs} and \code{sim} has to have the same length/dimension. \cr

Missing values in \code{obs} and/or \code{sim} can be removed before the computations, depending on the value of \code{na.rm}. \cr

Although \code{r} and \code{r2} have been widely used for model evaluation, these statistics are over-sensitive to outliers and insensitive to additive and proportional differences between model predictions and measured data (Legates and McCabe, 1999)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{me}}, \code{\link{mae}}, \code{\link{rmse}}, \code{\link{nrmse}}, \code{\link{pbias}},
\code{\link{pbiasfdc}}, \code{\link{rSD}}, \code{\link{NSeff}}, \code{\link{mNSeff}}, \code{\link{rNSeff}}, 
\code{\link{d}}, \code{\link{md}}, \code{\link{rd}}, \code{\link{cp}}, \code{\link{br2}}
}
\examples{
sim <- 1:10
obs <- 1:10
gof(sim, obs)

sim <- 2:11
obs <- 1:10
gof(sim, obs)

##################
# Loading daily streamflows of the Ega River (Spain), from 1961 to 1970
require(zoo)
data(EgaEnEstellaQts)
obs <- EgaEnEstellaQts

# Generating a simulated daily time series, initially equal to the observed series
sim <- obs 

# Getting the numeric goodness of fit for the "best" (unattainable) case
gof(sim=sim, obs=obs)

# Randomly changing the first 2000 elements of 'sim', by using a normal distribution with mean 10 and standard deviation equal to 1 (default of 'rnorm').
sim[1:2000] <- obs[1:2000] + rnorm(2000, mean=10)

# Getting the new numeric goodness of fit
gof(sim=sim, obs=obs)

# Getting the graphical representation of 'obs' and 'sim' along with the numeric goodness of fit 
\dontrun{
require(hydroTSM)
ggof(sim=sim, obs=obs)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }