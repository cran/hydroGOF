\name{ggof}
\Rdversion{1.1}
\alias{ggof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Graphical Goodness of Fit
}
\description{
Graphical comparison between two vectors (numeric, ts or zoo), with several numerical goodness of fit printed as a legend. \cr
Missing values in observed and/or simulated values can removed before the computations. \cr
}
\usage{
ggof(sim, obs, na.rm = TRUE, dates, date.fmt = "\%Y-\%m-\%d", 
     pt.style = "ts", ftype = "o",  FUN, 
     gof.leg = TRUE,  digits=2, 
     legend=c("Sim", "Obs"), leg.cex=1,
     tick.tstep = "months", lab.tstep = "years", lab.fmt,
     cal.ini=NA, val.ini=NA,
     main, xlab = "Time", ylab=c("Q, [m3/s]"),  
     col = c("blue", "black"), 
     cex = c(0.5, 0.5), cex.axis=1.2, cex.lab=1.2,
     lwd = c(1, 1), lty = c(1, 3), pch = c(1, 9), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{sim}{
numeric or zoo object with with simulated values
}
  \item{obs}{
numeric or zoo object with observed values
}
  \item{na.rm}{
a logical value indicating whether 'NA' should be stripped before the computation proceeds. \cr
When an 'NA' value is found at the i-th position in \code{obs} \bold{OR} \code{sim}, the i-th value of \code{obs} \bold{AND} \code{sim} are removed before the computation.
}
  \item{dates}{
character, factor or Date object indicating how to obtain the dates for the corresponding values in the \code{sim} and \code{obs} time series \cr
If \code{dates} is a character or factor, it is converted into Date class, using the date format specified by \code{date.fmt}
}
  \item{date.fmt}{
OPTIONAL. character indicating the format in which the dates are stored in \var{dates}, \code{cal.ini} and \code{val.ini}. See \code{format} in \code{\link[base]{as.Date}}. Default value is \kbd{\%Y-\%m-\%d} \cr
ONLY required when \code{class(dates)=="character"} or \code{class(dates)=="factor"} or when \code{cal.ini} and/or \code{val.ini} is provided.
}
  \item{pt.style}{
Character indicating if the 2 ts have to be plotted as lines or bars. When \code{ftype} is NOT \kbd{o}, it only applies for the annual values. Valid values are: \cr
-) \kbd{ts} : (default) each ts is plotted as a lines along the 'x' axis \cr
-) \kbd{bar}: the 2 series are plotted as a barplot. 
}
  \item{ftype}{
Character indicating how many plots are desired by the user. Valid values are: \cr
-) \kbd{o}  : only the original \code{sim} and \code{obs} time series are plotted \cr
-) \kbd{dm} : it assumes that \code{sim} and \code{obs} are daily time series and Daily and Monthly values are plotted  \cr
-) \kbd{ma} : it assumes that \code{sim} and \code{obs} are daily or monthly time series and Monthly and Annual values are plotted \cr
-) \kbd{dma}: it assumes that \code{sim} and \code{obs} are daily time series and Daily, Monthly and Annual values are plotted
}
  \item{FUN}{
OPTIONAL, ONLY required when \code{ftype} is in \code{c('dm', 'ma', 'dma')}. Function that have to be applied for transforming from daily to monthly or annual time step (e.g., for precipitation FUN MUST be \code{sum}, for temperature and flow time series, FUN MUST be \code{mean})
}
  \item{gof.leg}{
logical, indicating if several numerical goodness of fit have to be computed between \code{sim} and \code{obs}, and plotted as a legend on the graph. If \code{leg.gof=TRUE}, then \code{x} is considered as observed and \code{y} as simulated values (for some gof functions this is important).
}
  \item{digits}{
OPTIONAL, only used when \code{leg.gof=TRUE}. Numeric, representing the decimal places used for rounding the goodness-of-fit indexes.
}
  \item{legend}{
character of length 2 to appear in the legend.  
}
  \item{leg.cex}{
OPTIONAL. ONLY used when \code{leg.gof=TRUE}. Character expansion factor *relative* to current 'par("cex")'. Used for text, and provides the default for 'pt.cex' and 'title.cex'. Default value = 1
}
  \item{tick.tstep}{
character, indicating the time step that have to be used for putting the ticks on the time axis. Valid values are: \cr
-) \kbd{days},   \cr
-) \kbd{months}, \cr
-) \kbd{years}  
}
  \item{lab.tstep}{
character indicating the time step that have to be used for putting the labels on the time axis. Valid values are: \cr
-) \kbd{days},   \cr
-) \kbd{months}, \cr
-) \kbd{years}  
}
  \item{lab.fmt}{
Character indicating with the format to be used for the label of the axis. See \code{format} in \code{\link[base]{as.Date}}. If not specified, it will try \kbd{"\%Y-\%m-\%d"} when \code{lab.tstep=="days"}, \kbd{"\%b"} when \code{lab.tstep=="month"}, and \kbd{"\%Y"} when \code{lab.tstep=="year"}. 
}
  \item{cal.ini}{
OPTIONAL. Character, indicating the date in which the calibration period started. \cr
ONLY used for drawing a vertical red line at this date. 
}
  \item{val.ini}{
OPTIONAL. Character, the date in which the validation period started. \cr
ONLY used for drawing a vertical red line at this date. 
}
   \item{main}{
character representing the main title of the plot.
}
  \item{xlab}{
label for the 'x' axis.
}
  \item{ylab}{
label for the 'y' axis.
}
  \item{col}{
numeric, representing the colors of \code{sim} and \code{obs}
}
  \item{cex}{
numeric, representing the values controlling the size of text and symbols of 'x' and 'y' with respect to the default
}
  \item{cex.axis}{
numeric, representing the magnification to be used for the axis annotation relative to 'cex'. See \code{\link[graphics]{par}}.
}
  \item{cex.lab}{
numeric, representing the magnification to be used for x and y labels relative to the current setting of 'cex'. See \code{\link[graphics]{par}}.
}
  \item{lwd}{
vector with the line width of \code{sim} and \code{obs}
}
  \item{lty}{
numeric with the line type of \code{sim} and \code{obs}
}
  \item{pch}{
numeric with the type of symbol for \code{x} and \code{y}. (e.g., 1: white circle; 9: white rhombus with a cross inside)
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
Plots observed and simulated values in the same graph.

If \code{gof.leg=TRUE}, it computes the numerical values of: \cr
'me', 'mae', 'rms', 'nrms', 'PBIAS',  'RSR, 'rSD', 'NSeff', 'mNSeff', 'rNSeff', 'd', 'md, 'rd', 'cp', 'r', 'r.Spearman', 'R2', 'bR2'
}
\value{
\item{me}{Mean Error}
\item{mae}{Mean Absolute Error}
%%\item{mse}{Mean Squared Error}
\item{rmse}{Root Mean Square Error}
\item{nrmse}{Normalized Root Mean Square Error}
\item{PBIAS}{Percent Bias}
%%\item{pbfdc}{ PBIAS in the slope of the midsegment of the Flow Duration Curve}
\item{RSR}{Ratio of RMSE to the Standard Deviation of the Observations, RSR = rms / sd(obs). ( 0 <= RSR <= +Inf )}
\item{rSD}{Ratio of Standard Deviations, rSD = sd(sim) / sd(obs)}
\item{NSeff}{Nash-Sutcliffe Efficiency ( -Inf <= NSeff <= 1 )}
\item{mNSeff}{Modified Nash-Sutcliffe Efficiency}
\item{rNSeff}{Relative Nash-Sutcliffe Efficiency }
\item{d}{Index of Agreement ( 0 <= d <= 1 )}
\item{md}{Modified Index of Agreement}
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
%%\note{
%%
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gof}}, \code{\link{plot2}}
}
\examples{
obs <- 1:10
sim <- 2:11

\dontrun{
ggof(sim, obs)
}

##################
# Loading daily streamflows of the Ega River (Spain), from 1961 to 1970
require(zoo)
data(EgaEnEstellaQts)
obs <- EgaEnEstellaQts

# Generating a simulated daily time series, initially equal to the observed series
sim <- obs 

# Getting the numeric goodness of fit for the "best" (unattainable) case
gof(sim=sim, obs=obs)

# Randomly changing the first 2000 elements of 'sim', by using a normal distribution 
# with mean 10 and standard deviation equal to 1 (default of 'rnorm').
sim[1:2000] <- obs[1:2000] + rnorm(2000, mean=10)

# Getting the new numeric goodness of fit 
gof(sim=sim, obs=obs)

# Getting the graphical representation of 'obs' and 'sim' along with the numeric 
# goodness of fit for the daily and monthly time series 
\dontrun{
require(hydroTSM)
ggof(sim=sim, obs=obs, ftype="dm", FUN=mean)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line