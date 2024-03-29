# File rd.R
# Part of the hydroGOF R package, https://github.com/hzambran/hydroGOF ; 
#                                 https://cran.r-project.org/package=hydroGOF
#                                 http://www.rforge.net/hydroGOF/
# Copyright 2010-2024 Mauricio Zambrano-Bigiarini
# Distributed under GPL 2 or later

########################################
# 'rd': Relative Index of Agreement    #
########################################
# Started: April 15th, 2010            #
# Updates: 01-Jun-2011                 #
#          20-Jul-2022 ; 29-Jul-2022   #
#          20-Jan-2024                 #
########################################
# Ref
# 1) Krause, P., Boyle, D. P., and Base, F.: 
#    Comparison of different efficiency criteria for hydrological model assessment, Adv. Geosci., 5, 89-97, 2005 

# Relative Index of Agreement (Willmott et al., 1984) range from 0.0 to 1.0 
# and the closer to 1 the better the performance of the model 

# 'obs' : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim' : numeric 'data.frame', 'matrix' or 'vector' with simulated values

# 'Result': Modified Index of Agreement between 'sim' and 'obs'

# This index was developed to  be sensitive to systematic over- or 
# under-prediction, in particular during low flow conditions.

# This index quantify the difference between simulated and observed values 
# in a relative way, in order to significatively reduce the influence of 
# the absolute differences of high flows and to give more weight to 
# over- or under-prediction of low flows. At the same time, differences in
# low flows become more important, because they are looked in a relative way.

rd <-function(sim, obs, ...) UseMethod("rd")

rd.default <- function (sim, obs, na.rm=TRUE,
                        fun=NULL, ...,
                        epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                        epsilon.value=NA){ 

  if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
       is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
  ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")

  epsilon.type <- match.arg(epsilon.type)

  # index of those elements that are present both in 'x' and 'y' (NON- NA values)
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
     
    # Testing for zero values in obs, which leads to -Inf as result
    zero.index <- which(obs==0)
    if (length(zero.index > 0) ) 
      warning("'rd' can not be computed: some elements in 'obs' are zero !", call.=FALSE)
     
    # the next two lines are required for avoiding an strange behaviour 
    # of the difference function when sim and obs are time series.
    if ( !is.na(match(class(sim), c("ts", "zoo"))) ) sim <- as.numeric(sim)
    if ( !is.na(match(class(obs), c("ts", "zoo"))) ) obs <- as.numeric(obs)
     
    # Mean of the observed values
    Om <- mean(obs, na.rm=na.rm)
     
    denominator <- sum( ( ( abs(sim - Om) + abs(obs - Om) ) / Om )^2 )
     
    if ( (denominator != 0) & (!is.na(denominator)) ) {      
      rd <- 1 - ( sum( ( (obs - sim) / obs)^2 ) / denominator )     
    } else {
        rd <- NA
        warning("'sum( ( ( abs(sim-Om) + abs(obs-Om) ) / Om )^2 ) = 0', it is not possible to compute 'rd'")  
      } # ELSE end

  } else {
      rd <- NA
      warning("There are no pairs of 'sim' and 'obs' without missing values !")
    } # ELSE end
     
  return(rd) 
} # 'rd.default' end


rd.matrix <- function (sim, obs, na.rm=TRUE,
                       fun=NULL, ...,
                       epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                       epsilon.value=NA){ 

 # Checking that 'sim' and 'obs' have the same dimensions
 if ( all.equal(dim(sim), dim(obs)) != TRUE )
    stop( paste("Invalid argument: dim(sim) != dim(obs) ( [", 
          paste(dim(sim), collapse=" "), "] != [", 
          paste(dim(obs), collapse=" "), "] )", sep="") )

 rd <- rep(NA, ncol(obs))       
          
 rd <- sapply(1:ncol(obs), function(i,x,y) { 
                 rd[i] <- rd.default( x[,i], y[,i], na.rm=na.rm, fun=fun, ..., 
                                      epsilon.type=epsilon.type, epsilon.value=epsilon.value)
                 }, x=sim, y=obs )    
                     
  names(rd) <- colnames(obs)
  return(rd)
     
} # 'rd.matrix' end


rd.data.frame <- function (sim, obs, na.rm=TRUE,
                           fun=NULL, ...,
                           epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                           epsilon.value=NA){ 
 
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  rd.matrix(sim=sim, obs=obs, na.rm=na.rm, fun=fun, ..., 
            epsilon.type=epsilon.type, epsilon.value=epsilon.value)
     
} # 'rd.data.frame' end


################################################################################
# Author: Mauricio Zambrano-Bigiarini                                          #
################################################################################
# Started: 22-Mar-2013                                                         #
# Updates: 20-Jul-2022 ; 29-Jul-2022                                           #
################################################################################
rd.zoo <- function(sim, obs, na.rm=TRUE,
                   fun=NULL, ...,
                   epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                   epsilon.value=NA){ 
    
    sim <- zoo::coredata(sim)
    if (is.zoo(obs)) obs <- zoo::coredata(obs)
    
    if (is.matrix(sim) | is.data.frame(sim)) {
       rd.matrix(sim, obs, na.rm=na.rm, fun=fun, ..., 
                 epsilon.type=epsilon.type, epsilon.value=epsilon.value)
    } else NextMethod(sim, obs, na.rm=na.rm, fun=fun, ..., 
                      epsilon.type=epsilon.type, epsilon.value=epsilon.value)
     
  } # 'rd.zoo' end
