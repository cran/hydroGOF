# File br2.R
# Part of the hydroGOF R package, https://github.com/hzambran/hydroGOF
#                                 https://cran.r-project.org/package=hydroGOF
#                                 http://www.rforge.net/hydroGOF/ ;
# Copyright 2009-2024 Mauricio Zambrano-Bigiarini
# Distributed under GPL 2 or later

################################################################################
# 'br2': Weighted R2                                                           #
# Coef. of  determination multiplied by the coef. of the regression line       #
################################################################################
# Author: Mauricio Zambrano-Bigiarini                                          #
################################################################################
# Started: 27-Oct-2009                                                         #
# Updates: 11-Mar-2020                                                         #
#          16-Jan-2023 ; 29-Nov-2023                                           #
#          20-Jan-2024                                                         #
################################################################################

# This index allows accounting for the discrepancy in the magnitude of two signals
# under or overpredictions, (depicted by 'b') as well as their dynamics (depicted by R2).

# Krause, P., Boyle, D. P., and Base, F.: Comparison of different efficiency 
# criteria for hydrological model assessment, Adv. Geosci., 5, 89-97, 2005

# Krstic, G., Krstic, N. S., Zambrano-Bigiarini, M. (2016). 
# The br2-weighting Method for Estimating the Effects of Air Pollution on 
# Population Health. Journal of Modern Applied Statistical Methods, 15(2), 42. 
# doi:10.22237/jmasm/1478004000}

br2 <-function(sim, obs, ...) UseMethod("br2")
 
# 'obs'   : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim'   : numeric 'data.frame', 'matrix' or 'vector' with simulated values
# 'Result': weighted R2 between 'sim' and 'obs'

br2.default <- function(sim, obs, na.rm=TRUE, use.abs=FALSE, fun=NULL, ...,
                        epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                        epsilon.value=NA){

  if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
       is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
     ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")

  # the next two lines are required for avoiding an strange behaviour 
  # of the difference function when sim and obs are time series.
  if ( !is.na(match(class(sim), c("ts", "zoo"))) ) sim <- as.numeric(sim)
  if ( !is.na(match(class(obs), c("ts", "zoo"))) ) obs <- as.numeric(obs)

  # Checking 'epsilon.type'
  epsilon.type <- match.arg(epsilon.type)  

  # index of those elements that are present both in 'sim' and 'obs' (NON- NA values)
  vi <- valindex(sim, obs)
   
  if (length(vi) > 0) {	 
    # Filtering 'obs' and 'sim', selecting only those pairs of elements 
    # that are present both in 'x' and 'y' (NON- NA values)
    obs <- obs[vi]
    sim <- sim[vi]

    if (!is.null(fun)) {
      fun1 <- match.fun(fun)
      new  <- preproc(sim=sim, obs=obs, fun=fun1, ..., 
                      epsilon.type=epsilon.type, epsilon.value=epsilon.value)
      sim  <- new[["sim"]]
      obs  <- new[["obs"]]
    } # IF end     
     
    # Computing the linear regression between 'sim' and 'obs', 
    # forcing a zero intercept.
    x.lm <- lm(sim ~ obs - 1)
     
    # Getting the slope of the previous linear regression
    b <- as.numeric( coefficients(x.lm)["obs"]   )
     
    # computing the r2
    #r2 <- (.rPearson(sim, obs))^2 # this works only for linear models. 
    # https://github.com/hzambran/hydroGOF/issues/16#issue-1736556320
    r2 <- R2(sim=sim, obs=obs)

    if (!(use.abs)) {
      br2 <- ifelse(b <= 1, r2*abs(b), r2/abs(b))
    } else br2 <- ifelse(abs(b) <= 1, r2*abs(b), r2/abs(b))
   
  } else {
      br2 <- NA
      warning("There are no pairs of 'sim' and 'obs' without missing values !")
    } # ELSE end  
  
  return(br2)
     
} # 'br2' END
  

################################################################################
# 'br2': Weighted R2                                                           #
# Coef. of  determination multiplied by the coef. of the regression line       #
################################################################################
# Started: 27-Oct-2009                                                         #
# Updates: 11-Mar-2020                                                         #
#          16-Jan-2023                                                         #
################################################################################
br2.matrix <- function(sim, obs, na.rm=TRUE, use.abs=FALSE, fun=NULL, ...,
                       epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                       epsilon.value=NA){

    # Checking that 'sim' and 'obs' have the same dimensions
    if ( all.equal(dim(sim), dim(obs)) != TRUE )
      stop( paste("Invalid argument: dim(sim) != dim(obs) ( [", 
            paste(dim(sim), collapse=" "), "] != [", 
            paste(dim(obs), collapse=" "), "] )", sep="") )

    br2 <- rep(NA, ncol(obs))       
          
    br2 <- sapply(1:ncol(obs), function(i,x,y) { 
              br2[i] <- br2.default( x[,i], y[,i], na.rm=na.rm, use.abs=use.abs, 
                                     fun=fun, ..., 
                                     epsilon.type=epsilon.type,  
                                     epsilon.value=epsilon.value)
            }, x=sim, y=obs )            
           
    return(br2)
     
} # 'br2.matrix' END
  


################################################################################
# 'br2': Weighted R2                                                           #
# Coef. of  determination multiplied by the coef. of the regression line       #
################################################################################
# Started: 27-Oct-2009                                                         #
# Updates: 11-Mar-2020                                                         #
#          16-Jan-2023                                                         #
################################################################################  
br2.data.frame <- function(sim, obs, na.rm=TRUE, use.abs=FALSE, fun=NULL, ...,
                           epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                           epsilon.value=NA){

    sim <- as.matrix(sim)
    obs <- as.matrix(obs)

    br2.matrix(sim, obs, na.rm=na.rm, use.abs=use.abs, fun=fun, ..., 
               epsilon.type=epsilon.type, epsilon.value=epsilon.value)        
     
} # 'br2.data.frame' END
  
  
################################################################################
# Author: Mauricio Zambrano-Bigiarini                                          #
################################################################################
# Started: 22-Mar-2013                                                         #
# Updates: 11-Mar-2020                                                         #
################################################################################
br2.zoo <- function(sim, obs, na.rm=TRUE, use.abs=FALSE, fun=NULL, ...,
                    epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                    epsilon.value=NA){
    
    sim <- zoo::coredata(sim)
    if (is.zoo(obs)) obs <- zoo::coredata(obs)
    
    if (is.matrix(sim) | is.data.frame(sim)) {
       br2.matrix(sim, obs, na.rm=na.rm, use.abs=use.abs, fun=fun, ..., 
                  epsilon.type=epsilon.type, epsilon.value=epsilon.value)
    } else NextMethod(sim, obs, na.rm=na.rm, use.abs=use.abs, fun=fun, ..., 
                      epsilon.type=epsilon.type, epsilon.value=epsilon.value)
     
} # 'br2.zoo' end
