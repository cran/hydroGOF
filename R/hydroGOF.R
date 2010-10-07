######################################################################
#                'lib_GoodnessOfFit'                                 #
######################################################################
#  Library with several functions for assessing the goodness of fit  #
#  between observed and simulated values, the most commonly used in  #
#  Hydrological Modelling                                            #
#                                                                    #
#  Author      : Mauricio Zambrano Bigiarini                         #
#  Stared      : September 2008                                      #
#  Further dev.: Dec 2008, Jan, Feb, March, Sep, Oct 2009            # 
#  Version     : 0.1.0 : 07-Sep-2009                                 #
#  Version     : 0.1.1 : 06-Oct-2009                                 #
#  Version     : 0.1.2 : 29-Oct-2009                                 #
#  Version     : 0.1.3 : 01-Dec-2009                                 #
#  Version     : 0.2.0 : 07-Oct-2010                                 #
#  Version     : 0.2.1 : ongoing...                                  #
#  Last Update : 07-Jul-2010                                         #
######################################################################
# At the begining, only numerical vectors, matrix and data.frame were# 
# allowed as arguments, but on March 04th, 2009, the 'ts' and 'zoo'  #
# classes were added            

# Citations:
# 1) Boyle, D. P., H. V. Gupta, and S. Sorooshian (2000), Toward Improved Calibration of Hydrologic Models: Combining the Strengths of Manual and Automatic Methods, Water Resour. Res., 36(12), 3663–3674.                                      #
# 2) Krause, P., Boyle, D. P., and Bäse, F.: Comparison of different efficiency criteria for hydrological model assessment, Adv. Geosci., 5, 89-97, 2005.
# 3) Legates, D. R., and G. J. McCabe Jr. (1999), Evaluating the Use of "Goodness-of-Fit" Measures in Hydrologic and Hydroclimatic Model Validation, Water Resour. Res., 35(1), 233–241. 
# 4) Moriasi, D.N., Arnold, J.G., Van Liew, M.W., Bingner, R.L., Harmel, R.D., Veith, T.L. 2007. Model evaluation guidelines for systematic quantification of accuracy in watershed simulations. Transactions of the ASABE. 50(3):885-900. 
# 5) Kitanidis, P. K., and R. L. Bras (1980), Real-Time Forecasting With a Conceptual Hydrologic Model 2. Applications and Results, Water Resour. Res., 16(6), 1034–1044. 
# 6) J.E. Nash and J.V. Sutcliffe, River flow forecasting through conceptual models. Part 1: a discussion of principles, J. Hydrol. 10 (1970), pp. 282–290.
# 7) Yapo P. O., Gupta H. V., Sorooshian S., 1996. Automatic calibration of conceptual rainfall-runoff models: sensitivity to calibration data. Journal of Hydrology. v181 i1-4. 23-48


######################################################################## 
# External packages required:                                          #
# 1) 'zoo', by the 'ggof' function                                     #
# 2) 'hydroTSM', by the 'ggof' and 'plotbands'                         #
######################################################################## 

########################################################################
# Included functions:
# 01) .intersect    : elements that  belongs to 2 vectors
# 02) valindex : index of the elements that belongs to both vectors
# 03) ssq           : Sum of Squared Residuals
# 04) me            : Mean Error
# 05) mae           : Mean Absolute Error			
# 06) rmse          : Root Square Mean Error
# 07) nrmse         : Normalized Root Square Mean Error
# 08) rSD           : Ratio of Standard Deviations 
# 09) NSeff         : Nash-sutcliffe Efficiency
# 10) mNSeff        : Modified Nash-sutcliffe Efficiency (without the squares)
# 11) IoA           : Index of Agreement
# 12) PI            : Persistence Index
# 13) Pbias         : Percent Bias  
# 14) gof           : Several numerical performance indexes for comparing two vectors, matrix or data.frames
#                     It computes all the previous mentioned gof functions, 
#                     with the exception of: 'ssq', 'r.SD'
# 15) ggof          : Graphical performance comparison between two vectors (numeric, ts or zoo)   
# 16) plot2         : It plots 2 time series on the same graph. It is a wrapper for the 'plot.zoo' function  
# 17) br2           : weighted coef. of determination
# 18) mse           : Mean Squared Error 
# 19) rsr           : Ratio of RMSE to the Standard Deviation of the Observations   
# 20) plotbands     : It plots a ts of simulated values and two confidence bands, with optional plot of observations 
# 21) pfactor       : % of observations that are within the given uncertainty bounds
# 22) rfactor       : Average width of the given uncertainty bounds divided by the standard deviation of the observations   
# 23) pbiasfdc      : PBIAS in the slope of the midsegment of the Flow Duration Curve  
# 
########################################################################



###########################START of ex-lib_Plot.R#########################

############################################################
# 'plot2':     Plots 2 time series on the same graph       #
#              It is a wrapper for the 'plot.zoo' function #
############################################################
# Started on March 04, 2009  #
# May 2009                   #
##############################
 
# 'x', 'y'     : time series that will be plotted.
#                class(x) & class(y) must be 'ts' or 'zoo'
#                If gof.leg=TRUE, then 'x' is considered as simulated and 'y'
#                as observed values (for some gof functions this is important)
# 'plot.type'  : String that indicates if the 2 ts have to be ploted in the 
#                same window or in two different vertical ones
#                Valid values are:
#                -) "single"  : (default) superimposes the 2 ts on a single plot
#                -) "multiple": plots the 2 series on 2 multiple vertical plots 
# 'pt.style'   : String that indicates if the 2 ts have to be plotted as lines or bars
#                Valid values are:
#                -) "ts" : (default) each ts is ploted as a lines along the 'x' axis
#                -) "bar": the 2 series are plotted as a barplot. 
# 'var.names'  : string vector with the types (names) of variables being plotted, 
#                e.g, "Precipitation", "Temperature" or "Flow"
#                Only used for labelling the axes 
# 'var.units'  : string representing the measurement unit of the variable 
#                being plotted, e.g., "mm" for precipitation, "C" for temperature, 
#                and "m3/s" for flow. 
# 'tick.tstep': string indicating the time step that have to be used for 
#               putting the ticks ont he time axis. 
#               Possible values are: 'days', 'months', 'years' 
# 'lab.tstep' : string indicating the time step that have to be used for 
#               putting the labels ont he time axis. 
#               Possible values are: 'days', 'months', 'years' 
# 'col'       : vector with the colors of 'x' and 'y'
# 'lwd'       : vector with the line width of 'x' and 'y'
# 'lty'       : vector with the line type of 'x' and 'y'
# 'pch'       : vector with the type of symbol for 'x' and 'y'. 
#                1: whithe circle; 9: white rhombus with a cross inside
# 'cex'       : vector with the values controlling the size of text and 
#                symbols of 'x' and 'y' with respect to the default
# 'add'        : logical indicating if other plots will be added in further calls
#                to this function.
#                -) 'add=FALSE' => the plot and the legend are plotted on the same graph
#                -) 'add=TRUE'  => the legend is plotted in a new graph, usually
#                                  when called from another function (e.g.: 'ggof')
# 'xlab'       : label for the 'x' axis
# 'ylab'       : label for the 'y' axis 
# 'gof.leg'    : boolean indicating if several goodness of fit have to be 
#                computed between both ts, and ploted as legends on the graph.
#                If gof.leg=TRUE, then 'x' is considered as observed and 'y'
#                as simulated values (for some gof functions this is important)
# 'digits'     : OPTIONAL, only used when 'gof.leg=TRUE'. Decimal places used for rounding the goodness-of-fit indexes
# 'leg.cex'    : OPTIONAL. Used for the GoF legend. Character expansion factor *relative* to current
#                'par("cex")'.  Used for text, and provides the default 
#                for 'pt.cex' and 'title.cex'. Default value = 1
# cal.ini      : OPTIONAL. Character with the date in which the calibration period started.
#                ONLY used for drawing a vertical red line at this date. 
# val.ini      : OPTIONAL. Character with the date in which the validation period started.
#                ONLY used for drawing a vertical red line at this date. 
# date.fmt     : character indicating the format in which the dates entered are stored in 'cal.ini' adn 'val.ini'. Default value is "\%Y-\%m-\%d"
# 'cex.axis'   : magnification of axis annotation relative to 'cex'. See '?par'.
# 'cex.lab'    : Magnification to be used for x and y labels relative to the current setting of 'cex'. See '?par'.
                  
                
plot2 <- function (x, y, 
                   plot.type = "multiple", 
                   
                   tick.tstep= "months", 
                   lab.tstep= "years", 
                   
                   main, 
                   xlab="Time", 
                   ylab=c("x", "y"),
                   
                   cal.ini=NA, 
                   val.ini=NA, 
                   date.fmt="%Y-%m-%d",                   
                   
                   gof.leg = FALSE, 
                   gof.digits=2, 
                   
                   legend=ylab,
                   leg.cex=1,                       
                        
                   col = c("black","blue"),
                   
                   cex=c(0.5,0.5),
                   cex.axis=1.2,
                   cex.lab=1.2,
                   
                   lwd= c(1,1), 
                   lty= c(1,3), 
                   pch= c(1,9),   
                   
                   pt.style = "ts",
                   add=FALSE,                   
                   
                    ...) {
                   
  require(zoo)

  # requesting 'hydroTSM' package:'vector2zoo', 'drawxaxis'
  require(hydroTSM)

  # Checking that the user provided 'x'
  if ( missing(x) ) 
         stop("Missing argument: 'x'")
         
  # Checking that the user provided 'y'
  if ( missing(y) ) 
         stop("Missing argument: 'y'")
  
  # Checking that the user provided a valid argument for 'x'       
  if (is.na(match(class(x), c("integer", "numeric","ts", "zoo") ) ) ) 
         stop("Invalid argument: 'class(x)' must be in c('integer', 'numeric', 'ts', 'zoo')")
         
  # Checking that the user provided a valid argument for 'y'   
  if (is.na(match(class(y), c("integer", "numeric", "ts", "zoo") ) ) ) 
         stop("Invalid argument: 'class(y)' must be in c('integer', 'numeric', 'ts', 'zoo')")
         
  # Checking that the user provided a valid argument for 'plot.type'       
  if (is.na(match(plot.type, c("single", "multiple") ) ) ) 
         stop("Invalid argument: 'plot.type' must be in c('single', 'multiple')")
         
  # If the user wants to draw a legned, it checks that the type of plot is 'single'
  if (gof.leg & (plot.type == "multiple") )
    stop("Invalid argument: For drawing a legend, 'plot.type' must be 'single'")
         
  # Checking that the user provided a valid argument for 'pt.style'       
  if (is.na(match(pt.style, c("ts", "bar") ) ) ) 
         stop("Invalid argument: 'pt.style' must be in c('ts', 'bar')")
         
  # Checking that the user provided a valid argument for 'tick.tstep'       
  if (is.na(match(tick.tstep, c("days", "months", "years") ) ) ) 
         stop("Invalid argument: 'tick.tstep' must be in c('days', 'months', 'years')")
         
  # Checking that the user provided a valid argument for 'lab.tstep'       
  if (is.na(match(lab.tstep, c("days", "months", "years") ) ) ) 
         stop("Invalid argument: 'lab.tstep' must be in c('days', 'months', 'years')")
         
  # If 'x' is 'ts' or 'zoo' and 'y' is only a vector, y is transformed into 
  # the same class of 'x', with the same times
  if ( !is.na(match(class(x), c("ts", "zoo") ) ) & 
       !is.na(match(class(y), c("integer", "numeric") ) ) ) {  
  
    # class(time(x))== "Date" for 'daily' and 'monthly' time series
    # class(time(x))== "character" for 'annual' time series
    if ( class(time(x)) == "Date" ) {
        y <- vector2zoo(y, dates=time(x))
    } else if ( class(time(x)) == "character" ) {
        y <- vector2zoo(y, dates=time(x), date.fmt="%Y")
        time(x) <- time(y) #'annual' time series
    } # ELSE END
    
  } # IF END
         
  # Checking that the user provied the same length for 'sim' and 'obs'      
  #if ( length(x) != length(y) )  
  #       stop("Invalid argument: 'obs' and 'sim' must have the same length")
         
  # If the user didn't provide a title for the plot, the default is used 
  if ( missing(main) ) main <- "Observed vs Simulated"   
  
  # 'x' axis title: If the user didn't provide a title for the 'x' axis, a default string is used 
  if ( missing(xlab) ) { xlab <- "Time" }   
  
  # If the user provided a value for 'cal.ini', it is transformed into a Date class
  if ( !missing(cal.ini) ) 
    cal.ini <- as.Date(cal.ini, format=date.fmt)
    
  # If the user provided a value for 'val.ini', it is transformed into a Date class
  if ( !missing(val.ini) ) 
    val.ini <- as.Date(val.ini, format=date.fmt)

  
  # If the legend has to be plotted AND no other plots will be added
  # IF 'add' is TRUE, the layout of the screen is set up by the calling procedure (usually 'ggof')
  if (gof.leg & add==FALSE) {  
            def.par <- par(no.readonly = TRUE) # save default, for resetting...     
            #par(mar=c(5, 2, 4, 0.5) + 0.1)
            # Setting up the screen with 1 rows and 2 columns
            layout( matrix( c(1,1,1,1,1,1,1,1,1,2,2), ncol=11, byrow=TRUE) ) 
            par(mar=c(5, 4, 4, 0) + 0.1)    
            on.exit(par(def.par))      
  } # ELSE end    
  
  # If the legend will not be plotted, the marginns are set to 'almost' the default values
  if (!gof.leg) {  
        par(mar=c(5, 4, 4, 2) + 0.1) # default values are par(mar=c(5, 4, 4, 4) + 0.1)
  } # ELSE end    
  
  
  # If the type of plot is "time series"
  if (pt.style=="ts") {
  
    # If both time series have to be ploted int he same plot area
    if (plot.type == "single") {
    
      # Plotting the Observed Time Series
      # xaxt = "n": is for avoiding drawing the x axis
      # cbind(x, y) is a multivariate time series
      # 'screens = 1' can be used instead of 'plot.type="single"'
      # 'screens = c(1,2)' can be used for plotting each ts in a separate screen with the same time axis 
      plot.zoo( cbind(x, y), plot.type=plot.type, xaxt = "n", type=c("o","o"), 
               lwd=lwd, lty= lty, col= col, pch= pch, 
               cex = cex, cex.axis=cex.axis, cex.lab=cex.lab,
               main=main, xlab=xlab, ylab= ylab, ... )
               
      # If the user provided a value for 'cal.ini', a vertical line is drawn
      if ( !missing(cal.ini) ) {
        abline(v=cal.ini, col="red", lty=1, lwd=2)
      } # IF end
      
      # If the user provided a value for 'cal.ini', a vertical line is drawn
      if ( !missing(val.ini) ) {
        abline(v=val.ini, col="red", lty=1, lwd=2)
      } # IF end
               
      # Drawing a legend with 'Obs' vs 'Sim' 
      # y.intersp=0.5, is for the vertical spacin in the legend
      # bty="n" => no box around the legend
      # 'inset=0.03' is usefult when plot.type= "multiple" for having a nice margin to the legend
      legend("topright", legend=legend, y.intersp=0.8, inset=0.03,
             bty="n", cex = leg.cex, col = col, lwd= lwd, lty= lty, pch=pch )  
      
      # Drawing the 'x' axis
      # If the user provided, in some way, valid values for being used as dates, 
      # they will be used, if not, only a numeric index will be used
      if ( !is.na(match(class(x), c("ts", "zoo") ) ) | !is.na(match(class(y), c("ts", "zoo") ) ) ) {
  
        if ( !is.na(match(class(x), c("ts", "zoo") ) ) ) { 
          z <- x
        } else z <- y
    
        # Draws monthly ticks in the X axis, but labels only in years
        drawxaxis(z, tick.tstep=tick.tstep, lab.tstep= lab.tstep, cex.axis=cex.axis, cex.lab=cex.lab) 
    
      } else Axis(side = 1, labels = TRUE)
               
    } else  #plot.type == "multiple"  
          {       
            # all the work (mainly Time axis) is made automatically be the 'plot.zoo' function 
            plot.zoo( cbind(x, y), plot.type=plot.type, type=c("o","o"), 
                       lwd=lwd, lty= lty, col= col, pch= pch, 
                       cex = cex, cex.axis=cex.axis, cex.lab=cex.lab,
                       main=main, xlab=xlab, ylab= ylab,...)
                         
      } # ELSE end 
      
  } else if (pt.style=="bar") {
    
        # Creation of the table that will be plotted as barplot
        b <- rbind(coredata(x),coredata(y))
        
        # Giving the as name to each bar the YEAR, because the 
        # bar plot is thought for being used ONLY for annual time series
        colnames(b) <- format( time(x), "%Y")
        
        # Barplot  
        barplot(b, beside=TRUE, axis.lty=1, col=col, density=25, angle=c(45,-45), 
                main=main, xlab=xlab, ylab= ylab, legend.text=legend, 
                cex.axis=cex.axis, cex.lab=cex.lab, ...)
       
        # Index of the bar corresponding to 'cal.ini'.
        # It is necessary to multiply it by 3 because for each year there are 3 vertical lines
        # It is necessary to substract 2, for shifting the line form the 3 line to the first one
        cal.index <- 3*which(colnames(b) == format( cal.ini, "%Y")) - 2
        # If the user provided a value for 'cal.ini', a vertical line is drawn
        if ( !missing(cal.ini) ) {
         abline(v=cal.index, col="red", lty=1, lwd=2)
        } # IF end
       
        # Index of the bar corresponding to 'cal.ini'.
        # It is necessary to multiply it by 3 because for each year there are 3 vertical lines
        # It is necessary to substract 2, for shifting the line form the 3 line to the first one
        val.index <- 3*which(colnames(b) == format( val.ini, "%Y")) - 2
        # If the user provided a value for 'cal.ini', a vertical line is drawn
        if ( !missing(val.ini) ) {
          abline(v=val.index, col="red", lty=1, lwd=2)
        } # IF end
      
    }  # ELSE end        
  
  
  # If the Goodness of Fit indexes have to be computed and plotted:
  if (gof.leg & plot.type == "single" ) {
  
   gof.xy <- gof(sim=as.numeric(x), obs=as.numeric(y), do.spearman=FALSE, do.pbfdc=FALSE, digits=gof.digits, ...)
   
   legend.position <- "center"
   par( mar=c(0.5, 0.5, 0.5, 0.5) ) # mar=c(bottom, left, top, right). Default values are: mar=c(5,4,4,2) + 0.1)
   plot.new() 
          
   # 'inset':  The optional 'inset' argument specifies how far the legend is inset from the plot margins.  
   #           If a single value is given, it is used for both margins; 
   #           if two values are given, the first is used for 'x'-distance, the second for 'y'-distance.
	 
   legend(legend.position,  y.intersp=1.2, cex =leg.cex, # bty="n", #inset=0.01,   
          c( paste( "ME =", gof.xy["ME", 1], sep=" "),
             paste( "MAE =", gof.xy["MAE", 1], sep=" "),
             #paste( "MSE =", gof.xy["MSE", 1], sep=" "),
             paste( "RMSE =", gof.xy["RMSE", 1], sep=" "),
             paste( "NRMSE% =", gof.xy["NRMSE %", 1], sep=" "),
             paste( "PBIAS% =", gof.xy["PBIAS %", 1], sep=" "),
             #paste( "pbiasFDC% =", gof.xy["pbiasFDC %", 1], sep=" "),
             paste( "RSR =", gof.xy["RSR", 1], sep=" "),
             paste( "rSD =", gof.xy["rSD", 1], sep=" "),             
             paste( "NSeff =", gof.xy["NSeff", 1], sep=" "),
             paste( "mNSeff =", gof.xy["mNSeff", 1], sep=" "),
             paste( "rNSeff =", gof.xy["rNSeff", 1], sep=" "),
             paste( "d =", gof.xy["d", 1], sep=" "),
             paste( "md =", gof.xy["md", 1], sep=" "),
             paste( "rd =", gof.xy["rd", 1], sep=" "),
             #paste( "cp =", gof.xy["cp", 1], sep=" "),
             paste( "r =", gof.xy["r", 1], sep=" "),
             paste( "R2 =", gof.xy["R2", 1], sep=" "), 
             paste( "bR2 =", gof.xy["bR2", 1], sep=" ")             
            ), title="GoF's:", title.col="darkblue",
             bg="azure"
           )
         
  } #IF END
  
} # 'plot2' end


###########################END of ex-lib_Plot.R#########################

####################################################
# 'intersect': elements that  belongs to 2 vectors #
####################################################
#     19-Jan-2009   #
#####################
# 'x'     : vector (numerical or character)
# 'y'     : vector (numerical or character)
# 'Result': intersection between 'x' and 'y', i.e., only those elements that are
#           present both in 'x' and 'y'.
.intersect <- function(x, y) { 
 
 return( y[match(x, y, nomatch = 0)] )

} # '.intersect' END



#######################################################################
# 'valindex': index of the elements that belongs to both vectors #
#######################################################################
#     19-Jan-2009   #
#####################
# 'x'     : vector (numerical or character)
# 'y'     : vector (numerical or character)
# 'Result': index containing the position in 'x' and 'y' where both vectors 
#           have valid elements (NON- NA)
valindex <- function(obs, sim) {  

   if ( length(obs) != length(sim) ) 
	  stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !") 
	  
   valid.index.obs <- which( !is.na(obs) ) 
   valid.index.sim <- which( !is.na(sim) ) 
 
   return( .intersect(valid.index.obs, valid.index.sim) )
     
} # 'valindex.vec' END



###################################
# 'ssq': Sum of Squared Residuals #
###################################
#   04-Mar-2009; 06-Sep-09        #
###################################
# I don't think it is very useful, but I just added because it is required 
# for doing some SWAt comparisons

ssq <-function(sim, obs, ...) UseMethod("ssq")
 
# 'obs'   : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim'   : numeric 'data.frame', 'matrix' or 'vector' with simulated values
# 'Result': Sum of the Squared Residuals between 'sim' and 'obs', 
#           with squared measurement units of 'sim' and 'obs'
ssq.default <- function (sim, obs, na.rm=TRUE, ...){

     if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
          is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
         ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")    
      
     if ( length(obs) != length(sim) ) 
	    stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !") 
	 
	 ssq <- sum( (sim - obs)^2, na.rm= na.rm)   
     
     return(ssq)
     
  } # 'ssq' END
  

ssq.matrix <- function (sim, obs, na.rm=TRUE, ...){

    ssq <- colSums( (sim - obs)^2, na.rm = na.rm)          
           
    return(ssq)
     
  } # 'ssq.matrix' END
  
  
ssq.data.frame <- function (sim, obs, na.rm=TRUE, ...){

    sim <- as.matrix(sim)
	obs <- as.matrix(obs)
	
	ssq.matrix(sim, obs, na.rm=na.rm, ...)        
     
  } # 'ssq.data.frame' END
  
  

############################
# 'me': Mean Error         #
############################
#   15-Dic-2008; 06-Sep-09 #
############################
# 'obs'   : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim'   : numeric 'data.frame', 'matrix' or 'vector' with simulated values
# 'Result': Mean Error between 'sim' and 'obs', in the same units of 'sim' and 'obs' 

me <-function(sim, obs, ...) UseMethod("me")

me.default <- function (sim, obs, na.rm=TRUE, ...){

     
  if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
          is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
     ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")    
      
  if ( length(obs) != length(sim) ) 
	 stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !") 
		      
  me <- mean( sim - obs, na.rm = na.rm)           
     
  return(me)
     
  } # 'me.default' end
  
  
me.matrix <- function (sim, obs, na.rm=TRUE, ...){

   me <- colMeans( sim - obs, na.rm= na.rm)  
   
   return(me)
     
  } # 'me' end
  
  
me.data.frame <- function (sim, obs, na.rm=TRUE,...){

   sim <- as.matrix(sim)
   obs <- as.matrix(obs)
	
   me.matrix(sim, obs, na.rm=na.rm, ...)
     
  } # 'me.data.frame' end
  
  

##############################
# 'mae': Mean Absolute Error #
##############################
#   15-Dic-2008; 06-Sep-09   #
##############################
# 'obs'   : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim'   : numeric 'data.frame', 'matrix' or 'vector' with simulated values
# 'Result': Mean Absolute Error between 'sim' and 'obs', in the same units of 'sim' and 'obs' 

mae <-function(sim, obs, ...) UseMethod("mae")

mae.default <- function (sim, obs, na.rm=TRUE, ...){

  if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
          is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
     ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")    
      
  if ( length(obs) != length(sim) ) 
	 stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !")        
      
  mae <- mean( abs(sim - obs), na.rm = TRUE) 
               
  return(mae)
     
} # 'mae.default' end
  
  
mae.matrix <- function (sim, obs, na.rm=TRUE, ...){

  mae <- colMeans( abs(sim - obs), na.rm= na.rm)  
                 
  return(mae)
     
  } # 'mae.matrix' end
  
  
mae.data.frame <- function (sim, obs, na.rm=TRUE,...){

  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  mae.matrix(sim, obs, na.rm=na.rm, ...)  
     
} # 'mae.data.frame' end
  
  
##################################
# 'rmse': Root Mean Square Error #
##################################
#   15-Dic-2008; 06-Sep-09       #
##################################
# 'obs'   : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim'   : numeric 'data.frame', 'matrix' or 'vector' with simulated values
# 'Result': Root Mean Square Error between 'sim' and 'obs', in the same units of 'sim' and 'obs'

rmse <-function(sim, obs, ...) UseMethod("rmse")

rmse.default <- function (sim, obs, na.rm=TRUE, ...) {

  if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
          is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
     ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")    
      
  if ( length(obs) != length(sim) ) 
	 stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !")     
	 
  rmse <- sqrt( mean( (sim - obs)^2, na.rm = na.rm) )
           
  return(rmse)
     
} # 'rmse.default' end
  

rmse.matrix <- function (sim, obs, na.rm=TRUE, ...) {

   rmse <- sqrt( colMeans( (sim - obs)^2, na.rm = na.rm, ...) )          
           
   return(rmse)
     
} # 'rmse.matrix' end


rmse.data.frame <- function (sim, obs, na.rm=TRUE, ...) {

   rmse <- sqrt( colMeans( (sim - obs)^2, na.rm = na.rm, ...) )          
           
   return(rmse)
     
} # 'rmse.data.frame' end
  

##############################################
# 'nrmse': Normalized Root Mean Square Error #
############################################## 
#   15-Dic-2008; 06-Sep-09    #
############################### 
# 'obs'   : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim'   : numeric 'data.frame', 'matrix' or 'vector' with simulated values

# 'norm'  : character, indicating the value to be used to normalise the RMS. Valid values are:
#           -) 'sdobs' : standard deviation of observations.
#           -) 'maxmin': difference between maximum and minimum observed values

# 'Result': Normalized Root Mean Square Error between 'sim' and 'obs', 
#           when multiplied by 100 its units is %

nrmse <-function(sim, obs, ...) UseMethod("nrmse")
 
nrmse.default <- function (sim, obs, na.rm=TRUE, norm="sd", ...) {

    # Checking that the user provied a valid argument for 'norm'       
    if (is.na(match(norm, c("sd", "maxmin") ) ) ) 
       stop("Invalid argument: 'norm' must be in c('sd', 'maxmin')")
       
    if (norm=="sd") {
      cte <- sd(obs, na.rm=na.rm)
    } else if (norm=="maxmin") {
        cte <- ( max(obs, na.rm= na.rm) - min(obs, na.rm =na.rm) )
      } # ELSE end

     rmse <- rmse(sim, obs, na.rm) 
     
     if (max(obs, na.rm= na.rm) - min(obs, na.rm= na.rm) != 0) {
     
       nrmse <- rmse / cte
     
     } else stop("'obs' is constant, it is not possible to compute 'nrmse'")  
     
     return( round( 100*nrmse, 1) )
     
  } # 'nrmse.default' end
  
  
nrmse.matrix <- function (sim, obs, na.rm=TRUE, norm="sd", ...) {

  # Checking that the user provied a valid argument for 'norm'       
  if (is.na(match(norm, c("sdobs", "maxmin") ) ) ) 
     stop("Invalid argument: 'norm' must be in c('sd', 'maxmin')")

  nrmse <- rep(NA, ncol(obs))       
          
  nrmse <- sapply(1:ncol(obs), function(i,x,y) { 
                 nrmse[i] <- nrmse.default( x[,i], y[,i], na.rm=na.rm, norm=norm, ... )
               }, x=sim, y=obs )    
                     
  names(nrmse) <- colnames(obs)
  
  return(nrmse)
     
} # 'nrms.matrix' end


nrmse.data.frame <- function (sim, obs, na.rm=TRUE, norm="sd", ...) {

  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  nrmse.matrix(sim, obs, na.rm=na.rm, norm=norm, ...)
     
} # 'nrmse.data.frame' end
  
  
  
################################################
# 'rSD': Ratio of Standard Deviations          #
################################################
#   15-Dic-2008; 06-Sep-09    #
###############################
# SD(sim) / SD(obs)  
# 'obs'   : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim'   : numeric 'data.frame', 'matrix' or 'vector' with simulated values
# 'Result': Ratio of Standard Deviations  between 'sim' and 'obs', 
#           when multiplied by 100 its units is % 

rSD <-function(sim, obs, ...) UseMethod("rSD")

rSD.default <- function (sim, obs, na.rm=TRUE, ...){

     denominator <- sd(obs, na.rm = na.rm)
     
     if (denominator != 0) {
     
     rSD <- sd(sim, na.rm= na.rm) / sd(obs, na.rm= na.rm) 
     
     } else stop("'sd(obs)=0', it is not possible to compute 'rSD'")  
     
     return(rSD)
     
  } # 'rSD.default' end
  
  
rSD.matrix <- function (sim, obs, na.rm=TRUE, ...){ 
 
  rSD <- rep(NA, ncol(obs))       
          
  rSD <- sapply(1:ncol(obs), function(i,x,y) { 
                 rSD[i] <- rSD.default( x[,i], y[,i], na.rm=na.rm, ... )
               }, x=sim, y=obs )    
                     
  names(rSD) <- colnames(obs)
  
  return(rSD)
     
} # 'rSD.matrix' end


rSD.data.frame <- function (sim, obs, na.rm=TRUE, ...){ 
 
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  rSD.matrix(sim, obs, na.rm=na.rm, ...)
     
} # 'rSD.data.frame' end
  

########################################
# 'NSeff': Nash-sutcliffe Efficiency   #
########################################
# 15-Dic-2008   ; 06-Sep-09            #
########################################
# Nash-Sutcliffe efficiencies (Nash and Sutcliffe, 1970) range from -∞ to 1. 
# An efficiency of 1 (NSeff = 1) corresponds to a perfect match of modeled to the observed data. 
# An efficiency of 0 (NSeff = 0) indicates that the model predictions are as accurate
# as the mean of the observed data, whereas 
# an efficiency less than zero (-∞ < NSeff < 0) occurs when the observed mean is a better predictor than the model.
# Essentially, the closer the model efficiency is to 1, the more accurate the model is.  

# 'obs'   : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim'   : numeric 'data.frame', 'matrix' or 'vector' with simulated values
# 'Result': Nash-sutcliffe Efficiency between 'sim' and 'obs'

NSeff <-function(sim, obs, ...) UseMethod("NSeff")

NSeff.default <- function (sim, obs, na.rm=TRUE, ...){ 

   if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
          is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
     ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")      

   vi <- valindex(sim, obs)
     
   obs <- obs[vi]
   sim <- sim[vi]
     
   denominator <- sum( (obs - mean(obs))^2 )
     
   if (denominator != 0) {
      
     NS <- 1 - ( sum( (obs - sim)^2 ) / denominator )
     
   } else stop("'sum((obs - mean(obs))^2)=0' => it is not possible to compute 'NSeff'")  
     
   return(NS)
     
} # 'NSeff' end


NSeff.matrix <- function (sim, obs, na.rm=TRUE, ...){ 
 
  NS <- rep(NA, ncol(obs))       
          
  NS <- sapply(1:ncol(obs), function(i,x,y) { 
                 NS[i] <- NSeff.default( x[,i], y[,i], na.rm=na.rm, ... )
               }, x=sim, y=obs )    
                     
  names(NS) <- colnames(obs)
  
  return(NS)
     
} # 'NSeff.matrix' end


NSeff.data.frame <- function (sim, obs, na.rm=TRUE, ...){ 
 
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  NSeff.matrix(sim, obs, na.rm=na.rm, ...)
     
} # 'NSeff.data.frame' end


##################################################
# 'mNSeff': Modified Nash-sutcliffe Efficiency   #
##################################################
#    July 28th, 2009;  06-Sep-09                 #
##################################################
# Ref:
# Krause, P., Boyle, D. P., and Bäse, F.: Comparison of different efficiency criteria for hydrological model assessment, Adv. Geosci., 5, 89-97, 2005 
# Legates and McCabe, 1999. Evaluating the use of "goodness-of-fit" measures 
#                           in hydrologic and hydroclimatic model validation. 
#                           Water Resources Research. v35 i1. 233-241.

# Nash-Sutcliffe efficiency not "inflated" by squared values
# Essentially, the closer the model efficiency is to 1, the more accurate the model is.  

# 'obs' : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim' : numeric 'data.frame', 'matrix' or 'vector' with simulated values
# 'j'   : numeric, with the exponent to be used in the computation of the modified Nash-Sutcliffe effciency. The default value is j=1

# 'Result': Modified Nash-sutcliffe Efficiency between 'sim' and 'obs'

mNSeff <-function(sim, obs, ...) UseMethod("mNSeff")

mNSeff.default <- function (sim, obs, j=1, na.rm=TRUE, ...){ 

	 if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
          is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
     ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")   
     
     # Checking that the provided exponent is positive
     if (j < 0 ) stop("Invalid argument: 'j' must be positive")           
   
     vi <- valindex(sim, obs)
	 
	 obs <- obs[vi]
	 sim <- sim[vi]
	 
	 denominator <- sum( abs(obs - mean(obs))^j )
	 
	 if (denominator != 0) {
	  
	 NS1 <- 1 - ( sum( abs(obs - sim)^j ) / denominator )
	 
	 } else stop("'sum(abs(obs - mean(obs))^j)=0', it is not possible to compute 'mNSeff'")  
	 
	 return(NS1)
     
} # 'mNSeff.default' end


mNSeff.matrix <- function (sim, obs, j=1, na.rm=TRUE, ...){ 

  NS1 <- rep(NA, ncol(obs))       
          
  NS1 <- sapply(1:ncol(obs), function(i,x,y) { 
                 NS1[i] <- mNSeff.default( x[,i], y[,i], j, na.rm=na.rm, ... )
               }, x=sim, y=obs )    
                     
  names(NS1) <- colnames(obs)
  return(NS1)
     
} # 'mNSeff.matrix' end


mNSeff.data.frame <- function (sim, obs, j=1, na.rm=TRUE, ...){ 
 
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  mNSeff.matrix(sim, obs, j, na.rm=na.rm, ...)
     
} # 'mNSeff.data.frame' end



##################################################
# 'rNSeff': Relative Nash-sutcliffe Efficiency   #
##################################################
#    April 2010                                  #
##################################################
# Ref:
# Krause, P., Boyle, D. P., and Bäse, F.: Comparison of different efficiency criteria for hydrological model assessment, Adv. Geosci., 5, 89-97, 2005 
# Legates and McCabe, 1999. Evaluating the use of "goodness-of-fit" measures 
#                           in hydrologic and hydroclimatic model validation. 
#                           Water Resources Research. v35 i1. 233-241.

# Nash-Sutcliffe efficiency not "inflated" by squared values
# Essentially, the closer the model efficiency is to 1, the more accurate the model is.  

# 'obs' : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim' : numeric 'data.frame', 'matrix' or 'vector' with simulated values

# 'Result': Modified Nash-sutcliffe Efficiency between 'sim' and 'obs'

rNSeff <-function(sim, obs, ...) UseMethod("rNSeff")

rNSeff.default <- function (sim, obs, na.rm=TRUE, ...){ 

	 if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
          is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
     ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")      
   
     vi <- valindex(sim, obs)
	 
	 obs <- obs[vi]
	 sim <- sim[vi]
	 
	 denominator <- sum( ( ( obs - mean(obs) ) / mean(obs) )^2 )
	 
	 if (denominator != 0) {
	  
	 rNSeff <- 1 - ( sum( ( (obs - sim) / obs )^2 ) / denominator )
	 
	 } else stop("'sum( ( ( obs - mean(obs) ) / mean(obs) )^2 ) = 0', it is not possible to compute 'rNSeff'")  
	 
	 return(rNSeff)
     
} # 'rNSeff.default' end


rNSeff.matrix <- function (sim, obs, na.rm=TRUE, ...){ 

  rNSeff <- rep(NA, ncol(obs))       
          
  rNSeff <- sapply(1:ncol(obs), function(i,x,y) { 
                 rNSeff[i] <- rNSeff.default( x[,i], y[,i], na.rm=na.rm, ... )
               }, x=sim, y=obs )    
                     
  names(rNSeff) <- colnames(obs)
  return(rNSeff)
     
} # 'rNSeff.matrix' end


rNSeff.data.frame <- function (sim, obs, na.rm=TRUE, ...){ 
 
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  rNSeff.matrix(sim, obs, na.rm=na.rm, ...)
     
} # 'rNSeff.data.frame' end



########################################
# 'IoA': Index of Agreement            #
########################################
# December 18th, 2008;  06-Sep-09      #
########################################
# 1) Willmott, C.J. 1981. On the validation of models. Physical Geography, 2, 184-194
# 2) Willmott, C. J. (1984). "On the evaluation of model performance in physical geography." Spatial Statistics and Models, G. L. Gaile and C. J. Willmott, eds., 443-460.
# 3) Legates, D. R., and G. J. McCabe Jr. (1999), Evaluating the Use of "Goodness-of-Fit" Measures in Hydrologic and Hydroclimatic Model Validation, Water Resour. Res., 35(1), 233–241. 

# Index of Agreement (Willmott et al., 1984) range from 0.0 to 1.0 
# and the closer to 1 the better the performance of the model 

# 'obs'   : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim'   : numeric 'data.frame', 'matrix' or 'vector' with simulated values

# 'Result': Index of Agreement between 'sim' and 'obs'

d <-function(sim, obs, ...) UseMethod("d")

d.default <- function (sim, obs, na.rm=TRUE, ...){ 

     if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
          is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
     ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")

     # index of those elements that are present both in 'x' and 'y' (NON- NA values)
     vi <- valindex(sim, obs)
     
     # Filtering 'obs' and 'sim', selecting only those pairs of elements 
	 # that are present both in 'x' and 'y' (NON- NA values)
     obs <- obs[vi]
     sim <- sim[vi]
     
     # the next two lines are required for avoiding an strange behaviour 
     # of the difference function when sim and obs are time series.
     if ( !is.na(match(class(sim), c("ts", "zoo"))) ) sim <- as.numeric(sim)
     if ( !is.na(match(class(obs), c("ts", "zoo"))) ) obs <- as.numeric(obs)
     
     # Mean of the observed values
     Om <- mean(obs)
      
     denominator <- sum( ( abs(sim - Om) + abs(obs - Om)  )^2 )
     
     if (denominator != 0) {
      
       d <- 1 - ( sum( (obs - sim)^2 ) / denominator )
     
     } else stop("'sum((abs(sim-Om)+abs(obs-Om))^2)=0', it is not possible to compute 'IoA'")  
     
     return(d) 
     
} # 'd.default' end


d.matrix <- function (sim, obs, na.rm=TRUE, ...){ 

 d <- rep(NA, ncol(obs))       
          
 d <- sapply(1:ncol(obs), function(i,x,y) { 
                 d[i] <- d.default( x[,i], y[,i], na.rm=na.rm, ... )
                 }, x=sim, y=obs )    
                     
  names(d) <- colnames(obs)
  return(d)
     
} # 'd.matrix' end


d.data.frame <- function (sim, obs, na.rm=TRUE, ...){ 
 
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  d.matrix(sim, obs, na.rm=na.rm, ...)
     
} # 'd.data.frame' end



########################################
# 'md': Modified Index of Agreement    #
########################################
# April 07th, 2010                     #
########################################
# Ref
# 1) Krause, P., Boyle, D. P., and Bäse, F.: Comparison of different efficiency criteria for hydrological model assessment, Adv. Geosci., 5, 89-97, 2005 

# Index of Agreement (Willmott et al., 1984) range from 0.0 to 1.0 
# and the closer to 1 the better the performance of the model 

# 'obs' : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim' : numeric 'data.frame', 'matrix' or 'vector' with simulated values
# 'j'   : numeric, with the exponent to be used in the computation of the modified index  of agreement. The default value is j=1

# 'Result': Modified Index of Agreement between 'sim' and 'obs'

md <-function(sim, obs, ...) UseMethod("md")

md.default <- function (sim, obs, j=1, na.rm=TRUE, ...){ 

     if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
          is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
     ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")

     # Checking that the provided exponent is positive
     if (j < 0 ) stop("Invalid argument: 'j' must be positive")   

     # index of those elements that are present both in 'x' and 'y' (NON- NA values)
     vi <- valindex(sim, obs)
     
     # Filtering 'obs' and 'sim', selecting only those pairs of elements 
	 # that are present both in 'x' and 'y' (NON- NA values)
     obs <- obs[vi]
     sim <- sim[vi]
     
     # the next two lines are required for avoiding an strange behaviour 
     # of the difference function when sim and obs are time series.
     if ( !is.na(match(class(sim), c("ts", "zoo"))) ) sim <- as.numeric(sim)
     if ( !is.na(match(class(obs), c("ts", "zoo"))) ) obs <- as.numeric(obs)
     
     # Mean of the observed values
     Om <- mean(obs)
      
     denominator <- sum( ( abs(sim - Om) + abs(obs - Om)  )^j )
     
     if (denominator != 0) {
      
       d1 <- 1 - ( sum( ( abs(obs - sim) )^j ) / denominator )
     
     } else stop("'sum((abs(sim-Om)+abs(obs-Om))^j)=0', it is not possible to compute 'md'")  
     
     return(d1) 
     
} # 'md.default' end


md.matrix <- function (sim, obs, j=1, na.rm=TRUE, ...){ 

 d1 <- rep(NA, ncol(obs))       
          
 d1 <- sapply(1:ncol(obs), function(i,x,y) { 
                 d1[i] <- md.default( x[,i], y[,i], j, na.rm=na.rm, ... )
                 }, x=sim, y=obs )    
                     
  names(d1) <- colnames(obs)
  return(d1)
     
} # 'md.matrix' end


md.data.frame <- function (sim, obs, j=1, na.rm=TRUE, ...){ 
 
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  md.matrix(sim=sim, obs=obs, j=j, na.rm=na.rm, ...)
     
} # 'md.data.frame' end



########################################
# 'rd': Relative Index of Agreement    #
########################################
# April 15th, 2010                     #
########################################
# Ref
# 1) Krause, P., Boyle, D. P., and Bäse, F.: 
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

rd.default <- function (sim, obs, na.rm=TRUE, ...){ 

     if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
          is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
     ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")

     # index of those elements that are present both in 'x' and 'y' (NON- NA values)
     vi <- valindex(sim, obs)
     
     # Filtering 'obs' and 'sim', selecting only those pairs of elements 
	 # that are present both in 'x' and 'y' (NON- NA values)
     obs <- obs[vi]
     sim <- sim[vi]
     
     # the next two lines are required for avoiding an strange behaviour 
     # of the difference function when sim and obs are time series.
     if ( !is.na(match(class(sim), c("ts", "zoo"))) ) sim <- as.numeric(sim)
     if ( !is.na(match(class(obs), c("ts", "zoo"))) ) obs <- as.numeric(obs)
     
     # Mean of the observed values
     Om <- mean(obs)
      
     denominator <- sum( ( ( abs(sim - Om) + abs(obs - Om) ) / Om )^2 )
     
     if (denominator != 0) {
      
       rd <- 1 - ( sum( ( (obs - sim) / obs)^2 ) / denominator )
     
     } else stop("'sum( ( ( abs(sim-Om) + abs(obs-Om) ) / Om )^2 ) = 0', it is not possible to compute 'rd'")  
     
     return(rd) 
     
} # 'rd.default' end


rd.matrix <- function (sim, obs, na.rm=TRUE, ...){ 

 rd <- rep(NA, ncol(obs))       
          
 rd <- sapply(1:ncol(obs), function(i,x,y) { 
                 rd[i] <- rd.default( x[,i], y[,i], na.rm=na.rm, ... )
                 }, x=sim, y=obs )    
                     
  names(rd) <- colnames(obs)
  return(rd)
     
} # 'rd.matrix' end


rd.data.frame <- function (sim, obs, na.rm=TRUE, ...){ 
 
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  rd.matrix(sim=sim, obs=obs, na.rm=na.rm, ...)
     
} # 'rd.data.frame' end


 
########################################
# 'P': Coefficient of Persistence      #
########################################
# December 18th, 2008;  06-Sep-09      #
########################################
# Persistence Index (Kitadinis and Bras, 1980; Corradini et al., 1986) 
# is used to compare the model  performance agains a simple model using 
# the observed value of the previous day as the prediction for the current day.

#Kitanidis, P.K., and Bras, R.L. 1980. Real-time forecasting with a conceptual
#hydrologic model. 2. Applications and results. Water Resources Research,
#Vol. 16, No. 6, pp. 1034:1044.

# The coefficient of persistencec omparest he predictions of the model 
# with the predictions obtained by assuming that the process is a Wiener
# process(variance increasing linearly with time), in which case,
# the best estimate for the future is given by the latest measurement 
# (Kitadinis and Bras, 1980) 

# Persistence model efficiency (PME) is a normalized model evaluation statistic
# that quantifies the relative magnitude of the residual variance (noise)
# to the variance of the errors obtained by the use of a simple persistence 
# model 
# ("Ref: Moriasi, D.N., Arnold, J.G., Van Liew, M.W., Bingner, R.L., Harmel, 
#   R.D., Veith, T.L. 2007. Model evaluation guidelines for systematic 
#   quantification of accuracy in watershed simulations. 
#   Transactions of the ASABE. 50(3):885-900.. 

# PME ranges from 0 to 1, with PME = 1 being the optimal value. 
# PME values should be larger than 0.0 to indicate a minimally acceptable
# model performance (Gupta et al., 1999

# 'obs'   : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim'   : numeric 'data.frame', 'matrix' or 'vector' with simulated values
# 'Result': Persistence Index Efficiency between 'sim' and 'obs'

cp <-function(sim, obs, ...) UseMethod("cp")

cp.default <- function (sim, obs, na.rm=TRUE, ...){ 

     if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
          is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
     ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")

     # index of those elements that are present both in 'x' and 'y' (NON- NA values)
     vi <- valindex(sim, obs)
     
     # Filtering 'obs' and 'sim', selecting only those pairs of elements that are present both in 'x' and 'y' (NON- NA values)
     obs <- obs[vi]
     sim <- sim[vi]
     
     # the next two lines are required for avoiding an strange behaviour 
     # of the difference function when sim and obs are time series.
     if ( !is.na(match(class(sim), c("ts", "zoo"))) ) sim <- as.numeric(sim)
     if ( !is.na(match(class(obs), c("ts", "zoo"))) ) obs <- as.numeric(obs)
     
     # lenght of the data sets that will be ocnsidered for the ocmputations
     n <- length(obs)
      
     denominator <- sum( ( obs[2:n] - obs[1:(n-1)] )^2 )
     
     if (denominator != 0) {
      
     cp <- ( 1 - ( sum( (obs[2:n] - sim[2:n])^2 ) / denominator ) )
     
     } else stop("'sum((obs[2:n]-obs[1:(n-1))^2)=0', it is not possible to compute 'P'")  
     
     return(cp)
     
} # 'cp.default' end


cp.matrix <- function (sim, obs, na.rm=TRUE, ...){ 

  cp <- rep(NA, ncol(obs))       
          
  cp <- sapply(1:ncol(obs), function(i,x,y) { 
                 cp[i] <- cp.default( x[,i], y[,i], na.rm=na.rm, ... )
                 }, x=sim, y=obs )    
                     
   names(cp) <- colnames(obs)
     
   return(cp)
     
} # 'cp.matrix' end


cp.data.frame <- function (sim, obs, na.rm=TRUE, ...){ 
 
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  cp.matrix(sim, obs, na.rm=na.rm, ...)
     
} # 'cp.data.frame' end
  

##################################
# 'pbias': Percent Bias          #
##################################
#   03-Feb-2009;  06-Sep-09      #
##################################
# 'obs'   : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim'   : numeric 'data.frame', 'matrix' or 'vector' with simulated values
# 'Result': Percent Bias between 'sim' and 'obs', 
#           when multiplied by 100, its units is percentage
# Ref: Yapo P. O., Gupta H. V., Sorooshian S., 1996. 
#      Automatic calibration of conceptual rainfall-runoff models: 
#      sensitivity to calibration data. Journal of Hydrology. v181 i1-4. 23-48.

pbias <-function(sim, obs, ...) UseMethod("pbias")

pbias.default <- function (sim, obs, na.rm=TRUE, ...){

     if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
          is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
     ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")

     # index of those elements that are present both in 'x' and 'y' (NON- NA values)
     vi <- valindex(sim, obs)
     
     # Filtering 'obs' and 'sim', selecting only those pairs of elements 
	 # that are present both in 'x' and 'y' (NON- NA values)
     obs <- obs[vi]
     sim <- sim[vi]
     
     # lenght of the data sets that will be ocnsidered for the ocmputations
     n <- length(obs)
      
     denominator <- sum( obs )
     
     if (denominator != 0) {
      
       pbias <- 100 * ( sum( sim - obs ) / denominator )
     
     } else stop("'sum((obs)=0', it is not possible to compute 'pbias'")  
     
     return( round(pbias, 1) )
     
} # 'pbias.default' end
  
  
pbias.matrix <- function (sim, obs, na.rm=TRUE, ...){

   pbias <- rep(NA, ncol(obs))       
          
   pbias <- sapply(1:ncol(obs), function(i,x,y) { 
                 pbias[i] <- pbias.default( x[,i], y[,i], na.rm=na.rm, ... )
            }, x=sim, y=obs )        
                    
   return(pbias)
     
  } # 'pbias.matrix' end
  

pbias.data.frame <- function (sim, obs, na.rm=TRUE, ...){ 
 
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  pbias.matrix(sim, obs, na.rm=na.rm, ...)
     
} # 'pbias.data.frame' end  


##################################
#   rPearson;  27-Oct-2009       #
##################################
# Before Oct 27th 2009, this function was included in 'gof' function

# The 'r.Pearson' coefficient ranges from −1 to 1. 
# A value of 1 shows that a linear equation describes the relationship 
# perfectly and positively, with all data points lying on the same line 
# and with Y increasing with X. 
# A score of −1 shows that all data points lie on a single line but 
# that Y increases as X decreases. 
# A value of 0 shows that a linear model is not needed – that there 
# is no linear relationship between the variables.

.rPearson <-function(sim, obs, ...) UseMethod(".rPearson")

.rPearson.default <- function(sim, obs,...) {

  if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
          is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
     ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")
  
  rPearson <- cor(sim, obs, method="pearson", use="pairwise.complete.obs")      
  # if 'sim' and 'obs' were matrixs or data.frame, then the correlation
  # between observed and simulated values for each variable is given by the diagonal of 'r.Pearson' 
  
  #if ( is.matrix(r.Pearson) | is.data.frame(r.Pearson) ) {
  #r.Pearson        <- diag(r.Pearson)
  #}
  
  return(rPearson)
  
} # '.rPearson.default' end

.rPearson.matrix <- function (sim, obs, na.rm=TRUE, ...){

    rPearson <- rep(NA, ncol(obs))       
          
    rPearson <- sapply(1:ncol(obs), function(i,x,y) { 
                 rPearson[i] <- .rPearson.default( x[,i], y[,i], na.rm=na.rm, ... )
            }, x=sim, y=obs )            
           
    return(rPearson)
     
  } # '.rPearson.matrix' END
  
  
.rPearson.data.frame <- function (sim, obs, na.rm=TRUE, ...){

    sim <- as.matrix(sim)
	obs <- as.matrix(obs)
	
	.rPearson.matrix(sim, obs, na.rm=na.rm, ...)        
     
  } # '.rPearson.data.frame' END


#########################################################################
# 'br2': Weighted R2                                                    #
# Coef. of  determination multiplied by the coef. of the regression line#
#########################################################################
#   27-Oct-2009                                                         #
#########################################################################

# This index allows accounting for the discrepancy in the magnitude of two signals
# under or overpredictions, (depicted by 'b') as well as their dynamics (depicted by R2).

# Krause, P., Boyle, D. P., and Bäse, F.: Comparison of different efficiency 
# criteria for hydrological model assessment, Adv. Geosci., 5, 89-97, 2005

br2 <-function(sim, obs, ...) UseMethod("br2")
 
# 'obs'   : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim'   : numeric 'data.frame', 'matrix' or 'vector' with simulated values
# 'Result': weighted R2 between 'sim' and 'obs'

br2.default <- function (sim, obs, na.rm=TRUE, ...){

     if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
          is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
     ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")

     # index of those elements that are present both in 'x' and 'y' (NON- NA values)
     vi <- valindex(sim, obs)
     
     # Filtering 'obs' and 'sim', selecting only those pairs of elements 
     # that are present both in 'x' and 'y' (NON- NA values)
     obs <- obs[vi]
     sim <- sim[vi]
     
     # the next two lines are required for avoiding an strange behaviour 
     # of the difference function when sim and obs are time series.
     if ( !is.na(match(class(sim), c("ts", "zoo"))) ) sim <- as.numeric(sim)
     if ( !is.na(match(class(obs), c("ts", "zoo"))) ) obs <- as.numeric(obs)
     
     # Computing the linear regression between 'sim' and 'obs', 
     # forcing a zero intercept.
     x.lm <- lm(sim ~ obs - 1)
     
     # Getting the slope of the previous linear regression
     b <- as.numeric( coefficients(x.lm)["obs"]   )
     
     # computing the r2
     r2 <- (.rPearson(sim, obs))^2
     
     br2 <- ifelse(b <= 1, r2*abs(b), r2/abs(b))
     
     return(br2)
     
  } # 'br2' END
  

br2.matrix <- function (sim, obs, na.rm=TRUE, ...){

    br2 <- rep(NA, ncol(obs))       
          
    br2 <- sapply(1:ncol(obs), function(i,x,y) { 
                 br2[i] <- br2.default( x[,i], y[,i], na.rm=na.rm, ... )
            }, x=sim, y=obs )            
           
    return(br2)
     
  } # 'br2.matrix' END
  
  
br2.data.frame <- function (sim, obs, na.rm=TRUE, ...){

    sim <- as.matrix(sim)
    obs <- as.matrix(obs)

    br2.matrix(sim, obs, na.rm=na.rm, ...)        
     
  } # 'br2.data.frame' END
  
  
  
  
##################################
# 'mse': Mean Squared Error      #
##################################
#   03-Feb-2009;  06-Sep-09      #
##################################
# 'obs'   : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim'   : numeric 'data.frame', 'matrix' or 'vector' with simulated values
# 'Result': Mean Squared Error between 'sim' and 'obs'
# Ref: Yapo P. O., Gupta H. V., Sorooshian S., 1996. 
#      Automatic calibration of conceptual rainfall-runoff models: 
#      sensitivity to calibration data. Journal of Hydrology. v181 i1-4. 23-48.

mse <-function(sim, obs, ...) UseMethod("mse")

mse.default <- function (sim, obs, na.rm=TRUE, ...){

     if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
          is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
     ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")

     # index of those elements that are present both in 'x' and 'y' (NON- NA values)
     vi <- valindex(sim, obs)
     
     # Filtering 'obs' and 'sim', selecting only those pairs of elements 
     # that are present both in 'x' and 'y' (NON- NA values)
     obs <- obs[vi]
     sim <- sim[vi]
      
     mse <- mean( (sim - obs)^2, na.rm = na.rm)
     
     return( mse )
     
} # 'mse.default' end
  
  
mse.matrix <- function (sim, obs, na.rm=TRUE, ...){

   mse <- colMeans( (sim - obs)^2, na.rm= na.rm)        
                    
   return(mse)
     
  } # 'mse.matrix' end
  

mse.data.frame <- function (sim, obs, na.rm=TRUE, ...){ 
 
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  mse.matrix(sim, obs)
     
} # 'mse.data.frame' end  


########################################################################
# 'rsr': Ratio of RMSE to the Standard Deviation of the Observations   #
########################################################################
#   03-Feb-2010                  #
##################################
# 'obs'   : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim'   : numeric 'data.frame', 'matrix' or 'vector' with simulated values
# 'Result': Ratio of RMSE to the Standard Deviation of the Observations
#           It varies from 0 (its optimal value), which means zero RMSE and 
#           therefore a perfect model simulation, to +Inf. The lower the RSR, 
#           the better the model performance. Moriasi+al2007 suggest that 
#           a good performance is obtained for RSR < 0.7

# Ref: Moriasi, D.N., Arnold, J.G., Van Liew, M.W., Bingner, R.L., Harmel, 
#      R.D., Veith, T.L. 2007. Model evaluation guidelines for systematic 
#      quantification of accuracy in watershed simulations. 
#      Transactions of the ASABE. 50(3):885-900.

rsr <-function(sim, obs, ...) UseMethod("rsr")

rsr.default <- function (sim, obs, na.rm=TRUE, ...){

     if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo"))) |
          is.na(match(class(obs), c("integer", "numeric", "ts", "zoo")))
     ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo')")

     # index of those elements that are present both in 'x' and 'y' (NON- NA values)
     vi <- valindex(sim, obs)
     
     # Filtering 'obs' and 'sim', selecting only those pairs of elements 
     # that are present both in 'x' and 'y' (NON- NA values)
     obs <- obs[vi]
     sim <- sim[vi]
 
     #Root mean squared error
     rmse    <- rmse(sim=sim, obs=obs, na.rm=na.rm, ...)
     
     #Standard deviation of the observations
     sd.obs <- sd(obs, na.rm=na.rm)
     
     if ( sd.obs > 0 ) {
     
       rsr <- rmse / sd.obs
     
     } else stop("'sd(obs)=0', it is not possible to compute 'RSR'")  
     
     return( rsr )
     
} # 'rsr.default' end
  
  
rsr.matrix <- function (sim, obs, na.rm=TRUE, ...){

    rsr <- rep(NA, ncol(obs))       
          
    rsr <- sapply(1:ncol(obs), function(i,x,y) { 
                 rsr[i] <- rsr.default( x[,i], y[,i], na.rm=na.rm, ... )
            }, x=sim, y=obs )            
           
    return(rsr)  
     
  } # 'rsr.matrix' end
  

rsr.data.frame <- function (sim, obs, na.rm=TRUE, ...){ 
 
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  rsr.matrix(sim, obs, na.rm=na.rm, ...)
     
} # 'rsr.data.frame' end 



#########################################################################
# 'gof': Numerical Goodness of Fit between 'sim and 'obs'               #
#        Several performance indexes for comparing two vectors, matrix  #
#        or data.frames                                                 #
#########################################################################
#   15-Dic-2008 -> 03 Feb 2009; ;  06-Sep-09    #
#################################################

# It computes:
# 'me'        : Mean Error
# 'mae'       : Mean Absolute Error
# 'rms'       : Root Mean Square Error
# 'nrms'      : Normalized Root Mean Square Error
# 'r'         : Pearson Correltation coefficient ( -1 <= r <= 1 )
# 'r.Spearman': Spearman Correltation coefficient ( -1 <= r <= 1 ) 
# 'R2'        : Coefficient of Determination ( 0 <= r2 <= 1 )
#               Gives the proportion of the variance of one variable that
#               that is predictable from the other variable
# 'rSD'       : Ratio of Standard Deviations, rSD = SD(sim) / SD(obs)
# 'RSR'       : Ratio of the RMSE to the standard deviation of the observations
# 'NSeff'     : Nash-Sutcliffe Efficiency ( -Inf <= NSeff <= 1 )
# 'mNSeff'    : Modified Nash-Sutcliffe Efficiency
# 'd'         : Index of Agreement( 0 <= d <= 1 )
# 'md'        : Modified Index of Agreement( 0 <= md <= 1 )
# 'PI'        : Persistence Index ( 0 <= PI <= 1 ) 
# 'PBIAS'     : Percent Bias ( -1 <= PBIAS <= 1 )
# 'bR2'       : weighted coefficient of determination

gof <-function(sim, obs, ...) UseMethod("gof")

gof.default <- function (sim, obs, na.rm=TRUE, do.spearman=FALSE, do.pbfdc=FALSE, digits=2,...){

     ME     <- me(sim, obs, na.rm=na.rm)
     MAE    <- mae(sim, obs, na.rm=na.rm)
     MSE    <- mse(sim, obs, na.rm=na.rm)
     RMSE   <- rmse(sim, obs, na.rm=na.rm) 
     NRMSE  <- nrmse(sim, obs, na.rm=na.rm)
     RSR    <- rsr(sim, obs, na.rm=na.rm, ...)
     rSD    <- rSD(sim, obs, na.rm=na.rm)     
     PBIAS  <- pbias(sim, obs, na.rm=na.rm, ...)
     NSeff  <- NSeff(sim, obs, na.rm=na.rm, ...)
     mNSeff <- mNSeff(sim, obs, na.rm=na.rm, ...)
     rNSeff <- rNSeff(sim, obs, na.rm=na.rm, ...)
     d      <- d(sim, obs, na.rm=na.rm, ...)
     md     <- md(sim, obs, na.rm=na.rm, ...)
     rd     <- rd(sim, obs, na.rm=na.rm, ...)
     cp     <- cp(sim, obs, na.rm=na.rm, ...)
     r      <- .rPearson(sim, obs)
     bR2    <- br2(sim, obs, na.rm=na.rm, ...)     
     
     # 'r2' is the Coefficient of Determination
     # The coefficient of determination, r2, is useful because it gives the proportion of
     # the variance (fluctuation) of one variable that is predictable from the other variable.
     # It is a measure that allows us to determine how certain one can be in making
     # predictions from a certain model/graph.
     # The coefficient of determination is the ratio of the explained variation to the total
     # variation.
     # The coefficient of determination is such that 0 <  r2 < 1,  and denotes the strength
     # of the linear association between x and y. 
     R2 <- r^2
      
     if (do.spearman) {
       r.Spearman <- cor(sim, obs, method="spearman", use="pairwise.complete.obs") 
     
       # if 'sim' and 'obs' were matrixs or data.frame, then the correlation
       # between observed and simulated values for each variable is given by the diagonal of 'r.Pearson' 
       if ( is.matrix(r.Spearman) | is.data.frame(r.Spearman) ) {
         r.Spearman        <- diag(r.Spearman)
        } # IF end
        
     } # IF end
     
     if (do.pbfdc) { pbfdc  <- pbiasfdc(sim, obs, na.rm=na.rm, plot=FALSE, ...) }
     
     gof <- rbind(ME, MAE, MSE, RMSE, NRMSE, PBIAS, RSR, rSD, NSeff, mNSeff, rNSeff, d, md, rd, cp, r, R2, bR2)     
     
     rownames(gof)[5] <- "NRMSE %"
     rownames(gof)[6] <- "PBIAS %"    
     
     if (do.spearman) { gof <- rbind(gof, r.Spearman) }
     
     if (do.pbfdc) { 
       gof <- rbind(gof, pbfdc) 
       rownames(gof)[length(rownames(gof))] <- "pbiasFDC %"
     } # IF end
     
     # Rounding the final results, ofr avoiding scientific notation
     gof <- round(gof, digits)
     
     return(gof)
     
} # 'gof.default' end


gof.matrix <- function(sim, obs, na.rm=TRUE, do.spearman=FALSE, do.pbfdc=FALSE, digits=2, ...){

    # Temporal variable for some computations
    tmp <- gof(1:10,1:10)
    
    # Number of objective functions currently computed by gof
    ngof <- nrow(tmp) 
    
    # Name of the objective functions computed by 'gof'
    gofnames <- rownames(tmp)

    # Creating the matrix that will store the final results
    gof <- matrix(NA, ncol(obs), nrow=ngof)   
       
    # Computing the goodness-of-fit for each column of 'sim' and 'obs'      
    gof <- sapply(1:ncol(obs), function(i,x,y) { 
                 gof[, i] <- gof.default( x[,i], y[,i], na.rm=na.rm, do.spearman=do.spearman, do.pbfdc=FALSE, digits=digits, ... )
            }, x=sim, y=obs )            
     
    rownames(gof) <- gofnames
    colnames(gof) <- colnames(sim)
           
    return(gof)  
     
  } # 'gof.matrix' end
  

gof.data.frame <- function(sim, obs, na.rm=TRUE, do.spearman=FALSE, do.pbfdc=FALSE, digits=2,...){ 
 
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  gof.matrix(sim, obs, na.rm=na.rm, do.spearman=do.spearman, do.pbfdc=FALSE, digits=digits, ...)
     
} # 'gof.data.frame' end 
  
  

###########################################################################
# 'ggof': Graphical comparision bewtween two vectors (numeric, ts or zoo),#
#         with several numerical goodness of fit printed as a legend      #
###########################################################################
# External package required: 'zoo'                              #
# External library required: 'lib_TSM_in_Hydrological_modelling #
################################################################# 
#       Started:  03 Mar 2009           #
#  Apr, May 2009                        #
#########################################     


# 'var.types'	    : string representing the type of variable being plotted 
#                     (e.g., "Precipitation", "Temperature", "Flow",...). 
#                     ONLY used for labelling the axes 
# 'var.units'	    : string representing the measurement unit of the variable 
#                     being plotted (e.g., "mm" for precipitation, "C" for temperature,
#                     "m3/s" for flow,...).
#                     ONLY used for labelling the axes
# 'main'            : string representing the main title of the plot
# 'xlab'            : label for the 'x' axis
# 'ylab'            : label for the 'y' axis 
# 'ts.col'          : vector with the colors of 'sim' and 'obs'
# 'ts.lwd'          : vector with the line width of sim' and 'obs'
# 'ts.lty'          : vector with the line type of sim' and 'obs'
# 'ts.pch'          : vector with the type of symbol for 'x' and 'y'. 
#                     1: whithe circle; 9: white rhombus with a cross inside
# 'ts.cex'          : vector with the values controlling the size of text and 
#                     symbols of 'x' and 'y' with respect to the default
# 'ftype'           : string indicating the time frequency of the plots desired by the user. 
#                     Valid values are:
#                     -) 'o'  : only the original 'sim' and 'obs' time series are plotted
#                     -) 'dm' : it assumes that 'sim' and 'obs' are daily time series
#                               and Daily and Monthly values are plotted  
#                     -) 'ma' : it assumes that 'sim' and 'obs' are monthly time series
#                               and Monthly and Annual values are plotted
#                     -) 'dma': it assumes that 'sim' and 'obs' are daily time series
#                               and Daily, Monthly and Annual values are plotted 
# 'pt.style'        : String that indicates if the 2 ts have to be plotted as lines or bars
#                     When 'ftype' is NOT 'o', it only applies for the annual values
#                     Valid values are:
#                     -) "ts" : (default) each ts is ploted as a lines along the 'x' axis
#                     -) "bar": the 2 series are plotted as a barplot. 
# 'tick.tstep'      : string indicating the time step that have to be used for 
#                     putting the ticks ont he time axis. 
#                     Possible values are: 'days', 'months', 'years' 
# 'lab.tstep'       : string indicating the time step that have to be used for 
#                     putting the labels ont he time axis. 
#                     Possible values are: 'days', 'months', 'years'
# 'gof.leg'         : boolean indicating if several goodness of fit have to be 
#                     computed between both ts, and ploted as legends on the graph.
#                     If gof.leg=TRUE, then 'x' is considered as observed and 'y'
#                     as simulated values (for some gof functions this is important)
# 'digits'          : OPTIONAL, only used when 'gof.leg=TRUE'. Decimal places used for rounding the goodness-of-fit indexes
# 'leg.cex'         : Used for the GoF legend. Character expansion factor *relative* to current
#                     'par("cex")'.  Used for text, and provides the default 
#                     for 'pt.cex' and 'title.cex'. Default value = 0.7
# 'FUN'             : ONLY required when 'ftype' is in c('dm', 'ma', 'dma')
#                     Function that have to be applied for transforming from daily to monthly or annual time step
#                     For precipitation FUN MUST be "sum"
#                     For temperature and flow time series, FUN MUST be "mean"#             
# 'na.rm'           : Logical. ONLY matters when 'step.out' is "monthly' or 'annual'
#                     TRUE : the annual mean  value is computed considering only those values different from NA
#                     FALSE: if there is AT LEAST one NA within a year, the monthly mean value is NA 
# cal.ini           : OPTIONAL. Character with the date in which the calibration period started.
#                     ONLY used for drawing a vertical red line at this date. 
# val.ini           : OPTIONAL. Character with the date in which the validation period started.
#                     ONLY used for drawing a vertical red line at this date. 
# 'dates'           : Dates for the correponding values in the 'sim' and 'obs' time series
#                     If 'dates' is a factor, it have to be converted into 'Date' class, 
#                     using the date format  specified by 'date.fmt'
#                     If 'dates' is already of Date class, the number of dates
#                     must be equal to the number of elements in 'sim' and 'obs'
# date.fmt          : character indicating the format in which the dates entered are stored in 'cal.ini' adn 'val.ini'. Default value is "\%Y-\%m-\%d"
# 'cex.axis'        : magnification of axis annotation relative to 'cex'.
# 'cex.lab'         : Magnification to be used for x and y labels relative to the current setting of 'cex'. See '?par'.
                                               
      
ggof <- function (sim, obs, 
                  na.rm=TRUE, 
                  dates, 
                  date.fmt="%Y-%m-%d",

                  pt.style="ts",
                  ftype="o", 
                  FUN,
                  
                  gof.leg = TRUE, 
                  digits=2, 
                  
                  legend=c("Sim", "Obs"),
                  leg.cex=1,
                  
                  tick.tstep= "months", 
                  lab.tstep= "years",  
                  
                  cal.ini=NA, 
                  val.ini=NA,                
                  
                  main, 
                  xlab="Time", 
                  ylab=c("Q, [m3/s]"),  
                  
                  col= c("blue", "black"), 
                  
                  cex= c(0.5,0.5),
                  cex.axis=1.2,
                  cex.lab=1.2,
                  
                  lwd= c(1,1), 
                  lty= c(1,3), 
                  pch= c(1,9),                
                   
                   ...) {

  # requesting 'hydroTSM' package: 'sfreq', 'vector2zoo', 'daily2monthly', 'monthly2annual', 'daily2annual'
  require(hydroTSM)

  # Checking that the user provied a valid argument for 'sim'       
  if (is.na(match(class(sim), c("zoo", "numeric", "integer") ) ) ) 
         stop("Invalid argument: 'class(sim)' must be in c('zoo', 'numeric', 'integer')")
         
  # Checking that the user provied a valid argument for 'obs'       
  if (is.na(match(class(obs), c("zoo", "numeric", "integer") ) ) ) 
         stop("Invalid argument: 'class(obs)' must be in c('zoo', 'numeric', 'integer')")
         
  # Checking that the user provied the same length for 'sim' and 'obs'      
  if ( length(sim) != length(obs) )  
         stop(paste("Invalid argument: 'obs' and 'sim' must have the same length ! (", 
                   length(obs), "vs", length(sim), ")"  ,sep=" ") )
                   
  require(hydroTSM) # for using the 'sfreq' function
  # Checking that the user provied the same sampling frequency for 'sim' and 'obs',
  # when 'sim' and 'obs' are 'zoo' objects      
  if ( !is.na(match(class(obs), c("zoo") ) ) ) {
      if ( sfreq(sim) != sfreq(obs) ) {
         stop(paste("Invalid argument: 'obs' and 'sim' have different sampling frequency ! (", 
                   sfreq(obs), "vs", sfreq(sim), ")"  ,sep=" ") ) }
  } # IF end
         
  # Checking that the user provied a valid argument for 'ftype'       
  if (is.na(match(ftype, c("o", "dm", "ma", "dma") ) ) ) 
         stop("Invalid argument: 'ftype' must be in c('o', 'dm', 'ma, 'dma')")
         
  # Checking that the user provied a valid argument for FUN when 'ftype' involves monthly or annual values     
  if (!is.na(match(ftype, c("dm", "ma", "dma") ) ) & missing(FUN) ) 
         stop("Missing argument: 'FUN' must be provided when 'ftype' is in c('dm', 'ma, 'dma')")
         
  # If the user didn't provide a title for the plot, the default is used 
  if ( missing(main) ) main <- "Observations vs Simulations"
         
         
  # Requiring the Zoo Library (Zoo's ordered observations): 'is.zoo', 'as.zoo', and 'plot.zoo' functions
  require(zoo)
  
  # If the user provided values 'for 'dates'
  if (!missing(dates)) {
  
    # Checking that 'dates' have the same length than 'sim' ( and 'obs')      
    if ( length(dates) != length(sim) )  
         stop("Invalid argument: 'dates' and 'sim' must have the same length")
  
    # Checking that 'dates' have the right class
    if (is.na(match(class(dates), c("character", "factor", "Date")))) 
        stop("Invalid argument: 'class(dates)' must be in c('character', 'factor', 'Date')")
        
    # If 'dates' is a factor or character, it have to be converted into 'Date' class, 
    # using the date format  specified by 'date.fmt'
     if ( !is.na( match(class(dates), c("factor", "character") ) ) ) 
        dates <- as.Date(dates, format= date.fmt)    
    
    # If 'obs' is 'zoo' and the user provides the dates (probably new dates)
    if ( is.zoo(obs) ) { time(obs) <- dates }  
    # If 'sim' is 'zoo' and the user provides the dates  (probably new dates)
    if ( is.zoo(sim) ) { time(sim) <- dates }  
    
  } else if (!is.zoo(obs)) print("Note: You didn't provide dates, so only a numeric index will be used in the time axis.", quote=FALSE)      
 
  
  require(hydroTSM) # for using the 'vector2zoo' function 
  
  # If 'class(obs)' is not 'zoo' and the user provides the dates, then we turn it into a zoo class
  if ( !is.zoo(obs) & !missing(dates) ) { 
    obs <- vector2zoo(x=obs, dates=dates, date.fmt=date.fmt)        
  } # If 'class(obs)' is 'zoo' and 'dates' are missing, dates are extracted from 'obs'
    else if ( is.zoo(obs) & missing(dates) ) {  
      # class(time(x))== "Date" for 'daily' and 'monthly' time series
      # class(time(x))== "character" for 'annual' time series
      if ( class(time(obs)) == "Date" ) { dates <- time(obs) 
      } else if ( class(time(obs)) == "character" ) {  
             dates <- as.Date(time(obs), format="%Y") }      
    } #ELSE END
  
  # If 'class(sim)' is not 'zoo' and the user provides the dates, then we turn it into a zoo class
  if ( !is.zoo(sim) & !missing(dates) ) { 
    sim <- vector2zoo(x=sim, dates=dates, date.fmt=date.fmt) 
  # If 'class(sim)' is 'zoo' and 'dates' are missing, dates are extracted from 'sim'
  } else if ( is.zoo(sim) & is.zoo(obs) & missing(dates) ) {
      # class(time(x))== "Date" for 'daily' and 'monthly' time series
      # class(time(x))== "character" for 'annual' time series
      if ( class(time(sim)) == "Date" ) { dates <- time(obs) 
      } else if ( class(time(sim)) == "character" ) {  
             dates <- as.Date(time(sim), format="%Y") }
    } #ELSE END   
  
  
  #Plotting acoording to the 'ftype' value:  
  if (ftype == "o") {
        
   # Drawing the original time series against time
   plot2(x=sim, y=obs, plot.type="single",
         main= main, 
         col= col, lwd= lwd, lty=lty, pch=pch,
         xlab= xlab, ylab= ylab, pt.style= pt.style,
         add= FALSE,
         tick.tstep, lab.tstep, 
         cex = cex, cex.axis=cex.axis, cex.lab=cex.lab, 
         gof.leg = gof.leg, gof.digits=digits,
         legend=legend, leg.cex=leg.cex,
         cal.ini=cal.ini, val.ini=val.ini, date.fmt=date.fmt, ...)
         
  } else if (ftype=="dm") {
    
      if (sfreq(sim) != "daily") {      
        stop("Invalid argument: 'sim' has to have a 'daily' sampling frequency")       
      } else {
          # Generating a Monthly time series (Monthly mean of daily values):
          obs.monthly <- daily2monthly(obs, FUN, na.rm)
          sim.monthly <- daily2monthly(sim, FUN, na.rm)
          
          def.par <- par(no.readonly = TRUE) # save default, for resetting... 
          on.exit(par(def.par)) 
          
          # If the user wants a legend, the screen is splitted into 2 rows and 2 colums, 
          # where the proportion of width of the 1st column to the 2nd one is 9:2
          if (gof.leg) {           
            layout( matrix( c(1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,3,3,3,3,3,4,4), ncol=11, byrow=TRUE) ) 
          } else {
             # Setting up the screen with 2 rows and 1 column
             par(mfrow=c(2,1))
            } #ELSE end
          
          par(mar=c(5, 4, 4, 0) + 0.1) # mar=c(bottom, left, top, right). Default values are: mar=c(5,4,4,2) + 0.1)  
          # Drawing the original daily time series against time
          plot2(x=sim, y=obs, plot.type="single",
                main=paste("Daily", main, sep=" "), 
                tick.tstep=tick.tstep, lab.tstep= lab.tstep, 
                cex = cex, cex.axis=cex.axis, cex.lab=cex.lab, 
                col = col, lwd= lwd, lty=lty, pch=pch,  
                xlab= xlab, ylab= ylab, 
                pt.style= "ts", 
                add= TRUE,  
                gof.leg = gof.leg, gof.digits=digits,
                legend=legend, leg.cex=leg.cex,
                cal.ini=cal.ini, val.ini=val.ini, date.fmt=date.fmt, ... )
          
          # It is necessay to set up the margins again, after the previous call to plot2
          par(mar=c(5, 4, 4, 0) + 0.1) # mar=c(bottom, left, top, right). Default values are: mar=c(5,4,4,2) + 0.1)           
          # Drawing the Monthly time series against time
          plot2(x=sim.monthly, y=obs.monthly, plot.type="single",
                main=paste("Monthly", main, sep=" "), 
                tick.tstep=tick.tstep, lab.tstep= lab.tstep, 
                cex = cex, cex.axis=cex.axis, cex.lab=cex.lab, 
                col = col, lwd= lwd, lty=lty, pch=pch, 
                xlab= xlab, ylab= ylab, 
                pt.style= "ts", 
                add= TRUE, 
                gof.leg = gof.leg, gof.digits=digits,
                legend=legend, leg.cex=leg.cex,
                cal.ini=cal.ini, val.ini=val.ini, date.fmt=date.fmt, ... )
                   
            
        } # ELSE end
  } # ELSE if (ftype=="dm") END
  
  else if (ftype=="ma") {
  
    if  ( is.na( match( sfreq(sim), c("daily", "monthly") ) ) ) {      
      stop("Invalid argument: the sampling frequency of 'sim' has to be in c('daily', 'monthly'")       
    } else {
        if ( sfreq(sim) == "daily" ) {
           # Generating a Monthly time series (Monthly mean of daily values):
           obs <- daily2monthly(obs, FUN, na.rm)
           sim <- daily2monthly(sim, FUN, na.rm)
        } # IF end
        
        # Generating Annual time series (Annual mean of daily values)
        obs.annual <- monthly2annual(obs, FUN, na.rm, out.fmt="%Y-%m-%d")
        sim.annual <- monthly2annual(sim, FUN, na.rm, out.fmt="%Y-%m-%d")
        
        def.par <- par(no.readonly = TRUE) # save default, for resetting... 
        on.exit(par(def.par)) 
        
        # If the user wants a legend, the screen is splitted into 2 rows and 2 colums, 
        # where the proportion of width of the 1st column to the 2nd one is 9:2
        if (gof.leg) {     
          layout( matrix( c(1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,3,3,3,3,3,4,4), ncol=11, byrow=TRUE) )
        } else {
           # Setting up the screen with 2 rows and 1 column
           par(mfrow=c(2,1))
          } #ELSE end
        
        par(mar=c(5, 4, 4, 0) + 0.1) # mar=c(bottom, left, top, right). Default values are: mar=c(5,4,4,2) + 0.1)   
        # Drawing the Monthly time series against time
        plot2(x=sim, y=obs, plot.type="single",
              main=paste("Monthly", main, sep=" "), 
              tick.tstep=tick.tstep, lab.tstep= lab.tstep, 
              cex = cex, cex.axis=cex.axis, cex.lab=cex.lab, 
              col = col, lwd= lwd, lty=lty, pch=pch,
              xlab= xlab, ylab= ylab, pt.style= "ts", 
              add= TRUE, 
              gof.leg = gof.leg, gof.digits=digits,
              legend=legend, leg.cex=leg.cex,
              cal.ini=cal.ini, val.ini=val.ini, date.fmt=date.fmt, ... )
        
        # It is necessay to set up the margins again, after the previous call to plot2
        par(mar=c(5, 4, 4, 0) + 0.1)                
        # Drawing the Annual time series against time
        plot2(x=sim.annual, y=obs.annual, plot.type="single",
              main=paste("Annual", main, sep=" "), 
              tick.tstep="years", lab.tstep= "years", 
              cex = cex, cex.axis=cex.axis, cex.lab=cex.lab, 
              col = col, lwd= lwd, lty=lty, pch=pch, 
              xlab= xlab, ylab= ylab, pt.style= pt.style, 
              add= TRUE, 
              gof.leg = gof.leg, gof.digits=digits,
              legend=legend, leg.cex=leg.cex,
              cal.ini=cal.ini, val.ini=val.ini, date.fmt=date.fmt, ... )
      } # ELSE end      
   
  } # ELSE if (ftype=="ma") END
  
  else if (ftype=="dma") {
        
    if (sfreq(sim) != "daily") {      
      stop("Invalid argument: the 'sim' has to have a 'Daily' sampling frequency")  
           
    } else {       
          # Generating Monthly time series (Monthly mean of daily values):
          obs.monthly <- daily2monthly(obs, FUN, na.rm)
          sim.monthly <- daily2monthly(sim, FUN, na.rm)
          
          # Generating Annual time series (Annual mean of daily values)
          obs.annual <- daily2annual(obs, FUN, na.rm, out.fmt = "%Y-%m-%d")
          sim.annual <- daily2annual(sim, FUN, na.rm, out.fmt = "%Y-%m-%d")
          
          def.par <- par(no.readonly = TRUE) # save default, for resetting... 
          on.exit(par(def.par)) 
          
          # If the user wants a legend, the screen is splitted into 2 rows and 2 colums, 
          # where the proportion of width of the 1st column to the 2nd one is 9:2
          if (gof.leg) {   
            # Setting up a screen with 3 rows and 2 columns 
            layout( matrix( c(1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,3,3,3,3,3,4,4,5,5,5,5,5,5,5,5,5,6,6), ncol=11, byrow=TRUE) ) 
          } else {
             # Setting up the screen with 3 rows and 1 column
             par(mfrow=c(3,1))
            } #ELSE end  
          
          par(mar=c(5, 4, 4, 0) + 0.1) # mar=c(bottom, left, top, right). Default values are: mar=c(5,4,4,2) + 0.1) 
          # Drawing the original daily time series against time
          plot2(x=sim, y=obs, plot.type="single",
                main=paste("Daily", main, sep=" "), 
                tick.tstep=tick.tstep, lab.tstep= lab.tstep, 
                cex = cex, cex.axis=cex.axis, cex.lab=cex.lab, 
                col = col, lwd= lwd, lty=lty, pch=pch,
                xlab= xlab, ylab= ylab, pt.style= "ts", 
                add= TRUE, 
                gof.leg = gof.leg, gof.digits=digits,
                legend=legend, leg.cex=leg.cex,
                cal.ini=cal.ini, val.ini=val.ini, date.fmt=date.fmt, ... )
          
          # It is necessay to set up the margins again, after the previous call to plot2
          par(mar=c(5, 4, 4, 0) + 0.1) # mar=c(bottom, left, top, right). Default values are: mar=c(5,4,4,2) + 0.1)        
          # Drawing the Monthly time series against time
          plot2(x=sim.monthly, y=obs.monthly, plot.type="single",  
                main=paste("Monthly", main, sep=" "), 
                tick.tstep=tick.tstep, lab.tstep= lab.tstep, 
                cex = cex, cex.axis=cex.axis, cex.lab=cex.lab, 
                col = col, lwd= lwd, lty=lty, pch=pch, 
                xlab= xlab, ylab= ylab, pt.style= "ts", 
                add= TRUE, 
                gof.leg = gof.leg, gof.digits=digits,
                legend=legend, leg.cex=leg.cex,
                cal.ini=cal.ini, val.ini=val.ini, date.fmt=date.fmt, ... )
           
          # It is necessay to set up the margins again, after the previous call to plot2
          par(mar=c(5, 4, 4, 0) + 0.1) # mar=c(bottom, left, top, right). Default values are: mar=c(5,4,4,2) + 0.1)        
          # Drawing the Annual time series against time
          plot2(x=sim.annual, y=obs.annual, plot.type="single",
                  main=paste("Annual", main, sep=" "), 
                  tick.tstep="years", lab.tstep= "years", 
                  cex = cex, cex.axis=cex.axis, cex.lab=cex.lab, 
                  col = col, lwd= lwd, lty=lty, pch=pch,
                  xlab= xlab, ylab= ylab, pt.style= pt.style, 
                  add= TRUE, 
                  gof.leg = gof.leg, gof.digits=digits,
                  legend=legend, leg.cex=leg.cex,
                  cal.ini=cal.ini, val.ini=val.ini, date.fmt=date.fmt, ... )
      } # ELSE end
            
  } # ELSE if (ftype=="dma") END
  
} # 'ggof' end



########################################################################
# plotbands: Plot a ts with simulated values and two confidence bands  #
########################################################################
#	                   Date: 13-Oct-2009, 30-Jun-2010                  #
########################################################################
# 'x'        : ts or 'zoo' object with the observed values
# 'lband'    : ts or 'zoo' object with the values of the lower band
# 'uband'    : ts or 'zoo' object with the values of the upper band
# 'sim'      : OPTIONAL. ts or 'zoo' object with the simulated values for 'x'
# 'x.cex'    : character (or symbol) expansion for 'x' ts: a numerical vector.  This works as a multiple of par("cex")
# 'x.ty'     : character.  Indicates if the observed series have to be ploted as lines or points. Posiible values are:
#                -) "lines"  : the observed series are plotted as lines
#                -) "points" : the observed series are plotted as points
# 'lty'      : See '?plot.default'. The line type, see 'par'. 
# 'lwd'      : See '?plot.default'. The line width, see 'par'.}
# 'col'      : colors to be used for plotting the 'x' and 'sim' ts
# 'bands.col': color to be used for filling th area between the lower and upper band
# 'border'   : see '?polygon'. The color to draw the border of the uncertainty bands.  The default 'NA' omits the borders.
#             Use 'border' = 'NULL', to  use 'par("fg")'
# 'sim.col'  : OPTIONAL. color to be used for plotting the simulated ts
# 'main'     : an overall title for the plot: see 'title'
# 'xlab'     : a title for the x axis: see 'title'
# 'tick.tstep': string indicating the time step that have to be used for
#               putting the ticks ont he time axis.
#               Possible values are: 'days', 'months', 'years'
# 'lab.tstep' : string indicating the time step that have to be used for
#               putting the labels ont he time axis.
# 'cal.ini'   : OPTIONAL. Character with the date in which the calibration period started.
#               ONLY used for drawing a vertical red line at this date.
# 'val.ini'   : OPTIONAL. Character with the date in which the validation period started.
#               ONLY used for drawing a vertical red line at this date.
# 'date.fmt'  : character indicating the format in which the dates entered are stored in 'cal.ini' adn 'val.ini'. Default value is "\%Y-\%m-\%d"
# 'gof.leg'   : logical, indicating if the p-factor and r-factor have to be
#               computed and ploted as legends on the graph.
# 'gof.digits': OPTIONAL, numeric. Only used when 'gof.leg=TRUE'. 
#               Decimal places used for rounding p-factor and r-factor
# 'leg.cex'   : OPTIONAL. numeric. Used for the GoF legend. Character expansion factor *relative* to current
#               'par("cex")'.  Used for text, and provides the default
#               for 'pt.cex' and 'title.cex'. Default value = 1

# Example:
      
plotbands <- function(x, lband, uband, sim,
                      
                      dates,
                      date.fmt="%Y-%m-%d", 

                      gof.leg= TRUE, 
                      gof.digits=2,    
                      
                      legend=c("Obs", "Sim", "95PPU"),
                      leg.cex=1,
                      
                      bands.col="lightblue",
                      border= NA,

                      tick.tstep= "months", 
                      lab.tstep= "years",
                      
                      cal.ini=NA, 
                      val.ini=NA,                       
                      
                      
                      main="Confidence Bounds for 'x'",
                      xlab="Time",
                      ylab="Q, [m3/s]",
                      
                      col= c("black", "blue"),
                      type= c("lines", "lines"),

                      cex= c(0.5, 0.5),
                      cex.axis=1.2,
                      cex.lab=1.2,
                      
                      lwd=c(0.6, 1),
                      lty=c(3, 4),    
                      pch=c(1, 9),   
                      
                      ...) {
                    
    # requesting 'hydroTSM' package: 'sfreq', 'vector2zoo', 'drawxaxis'
    require(hydroTSM)

    # Checking  the class of 'x', 'lband', 'uband, and 'sim' (if provided)
    if ( is.na( match(class(x), c("zoo", "numeric", "integer") ) ) )
      stop("Invalid argument: 'class(x)' must be in c('zoo', 'numeric', 'integer')")
    if ( is.na( match(class(lband), c("zoo", "numeric", "integer") ) ) )
      stop("Invalid argument: 'class(lband)' must be in c('zoo', 'numeric', 'integer')")
    if ( is.na( match(class(uband), c("zoo", "numeric", "integer") ) ) )
      stop("Invalid argument: 'class(uband)' must be in c('zoo', 'numeric', 'integer')")      
    if ( !missing(sim) ) {
      if ( is.na( match(class(sim), c("zoo", "numeric", "integer") ) ) )
        stop("Invalid argument: 'class(sim)' must be in c('zoo', 'numeric', 'integer')")
    } # IF end    

    # Checking that the lenghts of 'lband' and 'uband' are equal between 
    # them, an equal to the lenght of 'x' and 'sim' (if provided)
    if ( length(lband) != length(uband) )
      stop("Invalid argument: 'length(lband)' is different from 'length(uband)'")
    if ( length(x) != length(uband) )
      stop("Invalid argument: 'length(x)' is different from the lengths of the bands")      
    if ( !missing(sim) ) {
      if ( length(x) != length(sim) )
      stop("Invalid argument: 'length(sim)' is different from 'length(x)'")
    } # IF end 
    
    # Length of the observed values and all the vectors provided
    L <- length(x)
      
      
    # Checking 'type'
    if ( length( which( !is.na( match((type), c("lines", "points") ) ) ) ) < 2)
      stop("Invalid argument: 'type' elements must be in c('lines', 'points')")
      
    # If the user provided a value for 'legend' here we verify that it has 3 elements
    if ( !missing(legend) )
      if ( length(legend) != 3 ) stop("Invalid argument: 'legend' must have 3 elements. e.g, c('obs', 'sim', 'PPU')")

    # If the user provided a value for 'cal.ini', it is transformed into a Date class
    if ( !missing(cal.ini) ) cal.ini <- as.Date(cal.ini, format=date.fmt)
    # If the user provided a value for 'val.ini', it is transformed into a Date class
    if ( !missing(val.ini) ) val.ini <- as.Date(val.ini, format=date.fmt)
      
    
    # Requiring the Zoo Library (Zoo's ordered observations): 'is.zoo', 'as.zoo', and 'plot.zoo' functions
    require(zoo)
    
    # If the user didn't provided the dates, but 'x' is a zoo object
    # dates are taken from 'x'
    if ( missing(dates) ) {
    
      if ( is.zoo(x) ) {
        # class(time(x))== "Date" for 'daily' and 'monthly' time series
        # class(time(x))== "character" for 'annual' time series
        if ( class(time(x)) == "Date" ) { dates <- time(x) 
        } else if ( class(time(x)) == "character" ) {  
             dates <- as.Date(time(x), format="%Y") 
          }  
      } else # If there is no way to obtain the dates
          print("Note: You didn't provide dates, so only a numeric index will be used in the time axis.", quote=FALSE)  
          
      # Checking that the dates of 'x', 'lband', 'uband' and 'sim' are equal ,
      # when they are zoo objects    
      if ( is.zoo(lband) & is.zoo(uband) ) 
        if  ( !all.equal( time(lband), time(uband) ) )
         stop("Invalid argument: time(lband) is different from time(uband)")       
      if ( is.zoo(x) & is.zoo(uband) ) 
        if  ( !all.equal( time(x), time(uband) ) )
          stop("Invalid argument: time(x) is different from the time of the bands")      
      if ( !missing(sim) ) {
        if ( is.zoo(x) & is.zoo(sim) ) 
          if  ( !all.equal( time(x), time(sim) ) )
            stop("Invalid argument: time(x) is different from the time of 'sim'")    
      } # IF end
          
    } # IF end
    
    # If the user provided 'dates', 
    # its length is checked against 'length(x)', and
    # the values of 'dates' are set to 'x', 'lband', 'uband' and 'sim' 
    # when they are zoo objects 
    if ( !missing(dates) )  { 
  
      # Checking that 'dates' have the same length than 'sim' ( and 'obs')      
      if ( length(dates) != length(x) )  
         stop("Invalid argument: 'dates' and 'x' must have the same length")
  
      # Checking that 'dates' have the right class
      if (is.na(match(class(dates), c("character", "factor", "Date")))) 
        stop("Invalid argument: 'class(dates)' must be in c('character', 'factor', 'Date')")
        
      # If 'dates' is a factor or character , it have to be converted into 'Date' class, 
      # using the date format  specified by 'date.fmt'
      if ( !is.na( match(class(dates), c("factor", "character") ) ) ) 
        dates <- as.Date(dates, format= date.fmt)   
    
      # If 'x', 'lband', 'uband' and 'sim' (when provided) are 'zoo' 
      # and the user provides 'dates' (probably new dates), 
      # the dates of the objects are changed to the new date
      if ( is.zoo(x) )     { time(x)     <- dates }  
      if ( is.zoo(lband) ) { time(lband) <- dates } 
      if ( is.zoo(uband) ) { time(uband) <- dates }  
      if ( !missing(sim) ) 
        if ( is.zoo(sim) ) { time(sim)   <- dates }  
        
      # If the class of 'x' 'lband', 'uband' and 'sim' (when provided) 
      # are not 'zoo' and the user provides the dates, 
      # then we turn them into a zoo objects
      if ( !is.zoo(x) )      x     <- vector2zoo(x=x, dates=dates, date.fmt=date.fmt) 
      if ( !is.zoo(lband) )  lband <- vector2zoo(x=lband, dates=dates, date.fmt=date.fmt) 
      if ( !is.zoo(uband) )  uband <- vector2zoo(x=uband, dates=dates, date.fmt=date.fmt) 
      if ( !missing(sim) )
        if ( !is.zoo(sim) )  sim   <- vector2zoo(x=sim, dates=dates, date.fmt=date.fmt)         
    
    }  # IF end
       
    
    # If 'x' is a zoo object, appropiate ticks and labels are set 
    # for the Time axis  
    if ( is.zoo(x) ) {    

        # If the user didn't provided a value for 'tick.tstep' and
        # the lenght of the daily ts is less than 1 year, the ticks in
        # the 'x' axis are placed by day; if larger than a year, they are placed by month
        if ( missing(tick.tstep) ) {    
          if ( sfreq(x)=="daily" ) {
            if ( (length(x) <= 366) ) {
              tick.tstep <- "days"
            } else tick.tstep <- "months"
          } else if ( sfreq(x)=="monthly" ) {
            if ( (length(x) <= 12) ) {
              tick.tstep <- "days"
            } else tick.tstep <- "months"
          } # ELSE end    
        } # IF end
    
        # If the user didn't provided a value for 'lab.tstep' and
        # the lenght of the daily ts is less than 1 year, the labels in
        # the 'x' axis are placed by month; if larger than a year, they are placed by year
        if ( missing(lab.tstep) ) {
          if ( sfreq(x)=="daily" ) {
            if (length(x) <=366) {
              lab.tstep <- "months"
            } else lab.tstep <- "years"
          } else if ( sfreq(x)=="monthly" ) {
            if (length(x) <=12) {
              lab.tstep <- "months"
            } else lab.tstep <- "years"
          }
        } # IF end
    
    } # IF end

    # Getting the position of the possible NA's
    na.index <- which(is.na(x))

    # Avoiding plotting the uncertainty bands for the Na's
    uband[na.index] <- uband[na.index-1]
    lband[na.index] <- lband[na.index+1]

    #uband[na.index] <- .5*( uband[na.index+1] + uband[na.index-1] )
    #lband[na.index] <- .5*( lband[na.index+1] + lband[na.index-1] )

    # Creating the 'x' values of the polygons of the bands
    if ( is.zoo(x) ) {
      t <- c( time(lband), rev(time(uband)) )
    } else t <- c( 1:L, L:1)

    # Creating the 'y' values of the polygons of the bands
    bands <- c(as.numeric(lband), rev(as.numeric(uband)) )
    
    # Computing 'ylim'
    if ( missing(sim) ) {
      ylim <- range(lband, uband, x, na.rm=TRUE)
    } else ylim <- range(lband, uband, x, sim, na.rm=TRUE)

    # Creating the plot, but without anything on it, for allowign the call to polygon
    plot(x, type="n", xaxt = "n", main=main, xlab=xlab, ylim=ylim, ...)

    # Plotting the polygons between the lower and upper bands
    polygon(t, bands, col=bands.col, border=border)

    # Draws monthly ticks in the X axis, but labels only in years
    if ( !is.na( match(class(x), c("zoo", "ts") ) ) ) {
      drawxaxis(x, tick.tstep=tick.tstep, lab.tstep=lab.tstep, cex.axis=cex.axis)
    } else Axis(side = 1, labels = TRUE)

    # Plotting the OBSERVED time series, over the polygons
    if (type[1] == "lines") {
      lines(x, cex= cex[1], col=col[1], lty=lty[1], lwd=lwd[1], pch=pch[1], ... )
    } else if (type[1] == "points") {
      points(x, cex= cex[1], col=col[1], lty=lty[1], lwd=lwd[1], pch=pch[1], ... )
      } # IF end
      
    # Plotting the SIMULATED time series, over the polygons
    if ( !missing(sim) ) {
        # Plotting the SIMULATED time series, over the polygons
        if (type[2] == "lines") {
          lines(sim, cex= cex[2], lty=lty[2], col=col[2], lwd=lwd[2], pch=pch[2], ... )
        } else if (type[2] == "points") {
          points(x, cex= cex[2], lty=lty[2], col=col[2], lwd=lwd[2], pch=pch[2], ... )
          } # IF end
    } # IF end

    # If the user provided a value for 'cal.ini', a vertical line is drawn
    if ( !missing(cal.ini) ) abline(v=cal.ini, col="red", lty=3, lwd=2)
    # If the user provided a value for 'val.ini', a vertical line is drawn
    if ( !missing(val.ini) ) abline(v=val.ini, col="red", lty=3, lwd=2)

    # Drawing a legend with 'Obs', 'Sim'  & 95PPU                        
    if ( missing(sim) ) {
      #legend <- c("Obs", "95PPU")
      legend <- c( legend[1], legend[3] ) 
      legend("topright", legend,  inset=0.03,
             bty="n", cex =0.9, col=c(col[1], bands.col), lwd=c(lwd[1], 0), lty=c(lty[1],0), pch=c(NA,15), pt.cex=3)
    } else {
      #legend <- c("Obs", "Sim", "95PPU") 
      legend("topright", legend, inset=0.03,
             bty="n", cex =0.9, col=c(col, bands.col), lwd=c(lwd, 0), lty=c(lty, 0), pch=c(NA,NA,15), pt.cex=3)
      } # ELSE end
    
    # Drawing a Legend with the p-factor and r-factor
    if (gof.leg) {
        legend("topleft",  y.intersp=1.2, cex =leg.cex, bty="n", #inset=0.01,
              c( paste( "p-factor =", round(pfactor(x, lband=lband, uband=uband), gof.digits), sep=" "),
                 paste( "r-factor =", round(rfactor(x, lband=lband, uband=uband), gof.digits), sep=" ") )
               )
    } # IF end

} # 'plotbands' END




##########################################################################
# pfactor: % of observations that are within the given uncertainty bounds#
##########################################################################
#	                       Date: 24-Jan-2010                             #
##########################################################################
# 'x'        : ts or 'zoo' object with the simulated values
# 'lband'    : ts or 'zoo' object with the values of the lower uncertainty bound
# 'uband'    : ts or 'zoo' object with the values of the upper uncertainty bound

# Ideally, i.e., with a combination of model structure and parameter values 
# that perfectly represents the catchment under study, and in absence of 
# measurement errors and other additional sources of uncertainty, all the 
# simulated values should be in a perfect match with the observations, 
# leading to a P-factor equal to 1, and an R-factor equal to zero. 
# However, in real-world applications we aim at encompassing as much 
# observations as possible within the given uncertainty bounds 
# (P-factor close to 1) while keeping the width of the uncertainty bounds 
# as small as possible (R-factor close to 0), in order to avoid obtaining 
# a good bracketing of observations at expense of uncertainty bounds too 
# wide to be informative for the decision-making process.

# Refs:
# Abbaspour, K. C., M. Faramarzi, S. S. Ghasemi, and H. Yang (2009), 
# Assessing the impact of climate change on water resources in Iran, 
# Water Resour. Res., 45(10), W10,434, doi:10.1029/2008WR007615. \cr

# Abbaspour, K. C., J. Yang, I. Maximov, R. Siber, K. Bogner, J. Mieleitner, 
# J. Zobrist, and R. Srinivasan (2007), Modelling hydrology and water quality 
# in the pre-alpine/alpine Thur watershed using SWAT, Journal of Hydrology, 
# 333(2-4), 413–430, doi:10.1016/j.jhydrol.2006.09.014.

# Schuol, J., K. Abbaspour, R. Srinivasan, and H. Yang (2008b), 
# Estimation of freshwater availability in the West African sub-continent 
# using the SWAT hydrologic model, Journal of Hydrology, 352(1-2), 30, 
# doi:10.1016/j.jhydrol.2007.12.025. \cr

# Abbaspour, C., Karim (2007), User manual for SWAT-CUP, SWAT calibration 
# and uncertainty analysis programs, 93pp, Eawag: Swiss Fed. Inst. of Aquat. Sci. and 
# Technol. Dubendorf, Switzerland, Available at http://www.eawag.ch/organisation/abteilungen/siam/software/swat/index_EN

pfactor <- function(x, ...) UseMethod("pfactor")

pfactor.default <- function(x, lband, uband, na.rm=TRUE, ...)  {

    # Just in case 'some' of them be ts or 'zoo'
    x     <- as.numeric(x)
    lband <- as.numeric(lband)
    uband <- as.numeric(uband)
        
    # Getting the row index in 'q95' of all the observations that are within L95PPU and U95PPU
    within.index <- which((lband <= x) & (x <= uband) )
     
    # Getting the best simulated streamflows (skipping days withoud measurements)
    pfactor <- length( within.index ) / length( x ) 
    
    return(pfactor)

} # 'pfactor.default' end


pfactor.matrix <- function (x, lband, uband, na.rm=TRUE, ...){

    pfactor <- rep(NA, ncol(x))       
          
    pfactor <- sapply(1:ncol(x), function(i,x,l,u) { 
                 pfactor[i] <- pfactor.default( x[,i], l[,i], u[i], na.rm=na.rm, ... )
            }, x=x, l=lband, u=uband )            
           
    return(pfactor)  
     
  } # 'pfactor.matrix' end
  

pfactor.data.frame <- function (x, lband, uband, na.rm=TRUE, ...){ 
 
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  pfactor.matrix(x, lband, uband, na.rm=na.rm, ...)
     
} # 'pfactor.data.frame' end 


#########################################################################
# rfactor: average width of the given uncertainty bounds divided by the #
#          standard deviation of the observations                       #
#########################################################################
#	                       Date: 24-Jan-2010                            #
#########################################################################
# 'x'        : ts or 'zoo' object with the simulated values
# 'lband'    : ts or 'zoo' object with the values of the lower band
# 'uband'    : ts or 'zoo' object with the values of the upper band

# Ideally, i.e., with a combination of model structure and parameter values 
# that perfectly represents the catchment under study, and in absence of 
# measurement errors and other additional sources of uncertainty, all the 
# simulated values should be in a perfect match with the observations, 
# leading to a P-factor equal to 1, and an R-factor equal to zero. 
# However, in real-world applications we aim at encompassing as much 
# observations as possible within the given uncertainty bounds 
# (P-factor close to 1) while keeping the width of the uncertainty bounds 
# as small as possible (R-factor close to 0), in order to avoid obtaining 
# a good bracketing of observations at expense of uncertainty bounds too 
# wide to be informative for the decision-making process.

# Refs:
# Abbaspour, K. C., M. Faramarzi, S. S. Ghasemi, and H. Yang (2009), 
# Assessing the impact of climate change on water resources in Iran, 
# Water Resour. Res., 45(10), W10,434, doi:10.1029/2008WR007615. \cr

# Abbaspour, K. C., J. Yang, I. Maximov, R. Siber, K. Bogner, J. Mieleitner, 
# J. Zobrist, and R. Srinivasan (2007), Modelling hydrology and water quality 
# in the pre-alpine/alpine Thur watershed using SWAT, Journal of Hydrology, 
# 333(2-4), 413–430, doi:10.1016/j.jhydrol.2006.09.014.

# Schuol, J., K. Abbaspour, R. Srinivasan, and H. Yang (2008b), 
# Estimation of freshwater availability in the West African sub-continent 
# using the SWAT hydrologic model, Journal of Hydrology, 352(1-2), 30, 
# doi:10.1016/j.jhydrol.2007.12.025. \cr

# Abbaspour, C., Karim (2007), User manual for SWAT-CUP, SWAT calibration 
# and uncertainty analysis programs, 93pp, Eawag: Swiss Fed. Inst. of Aquat. Sci. and 
# Technol. Dubendorf, Switzerland, Available at http://www.eawag.ch/organisation/abteilungen/siam/software/swat/index_EN

rfactor <- function(x, ...) UseMethod("rfactor")

rfactor.default <- function(x, lband, uband, na.rm=TRUE, ...)  {

    # Just in case 'some' of them be ts or 'zoo'
    x     <- as.numeric(x)
    lband <- as.numeric(lband)
    uband <- as.numeric(uband)
        
    # Getting the row index in 'q95' of all the observations that are within L95PPU and U95PPU
    dfactor <- mean(uband - lband, na.rm=na.rm)
     
    # Getting the best simulated streamflows (skipping days withoud measurements)
    rfactor <- dfactor / sd(x, na.rm=na.rm) 
    
    return( rfactor )

} # 'rfactor' end

rfactor.matrix <- function (x, lband, uband, na.rm=TRUE, ...){

    rfactor <- rep(NA, ncol(x))       
          
    rfactor <- sapply(1:ncol(x), function(i,x,l,u) { 
                 rfactor[i] <- rfactor.default( x[,i], l[,i], u[i], na.rm=na.rm, ... )
            }, x=x, l=lband, u=uband )            
           
    return(rfactor)  
     
  } # 'rfactor.matrix' end
  

rfactor.data.frame <- function (x, lband, uband, na.rm=TRUE, ...){ 
 
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  rfactor.matrix(x, lband, uband, na.rm=na.rm, ...)
     
} # 'rfactor.data.frame' end 


##############################################################################
# 'pbiasfdc': PBIAS in the slope of the midsegment of the Flow Duration Curve #
##############################################################################
#   03-Feb-2010                  #
##################################
# 'obs'   : numeric 'data.frame', 'matrix' or 'vector' with observed values
# 'sim'   : numeric 'data.frame', 'matrix' or 'vector' with simulated values
# 'lQ.thr': numeric, used to classify low flows. All the streamflows with a probaility of exceedence larger or equal to 'lQ.thr' are classified as low flows
# 'hQ.thr': numeric, used to classify high flows. All the streamflows with a probaility of exceedence lower or equal to 'hQ.thr' are classified as high flows
# 'plot'  : a logical value indicating if the flow duration curves corresponding to 'obs' and 'sim' have to be  plotted or not.

# 'Result': Percent Bias in the slope of the midsegment of the flow duration curve [%]
#           It measures the vertical soil moisture redistribution

# Ref:  Yilmaz, K. K., H. V. Gupta, and T. Wagener  (2008), 
#       A process-based diagnostic approach to model evaluation: 
#       Application to the NWS distributed hydrologic model, 
#       Water Resour. Res., 44, W09417, doi:10.1029/2007WR006716.


pbiasfdc <-function(sim, obs, ...) UseMethod("pbiasfdc")

pbiasfdc.default <- function (sim, obs, lQ.thr=0.7, hQ.thr=0.2, na.rm=TRUE, plot=TRUE, verbose=FALSE, ...){

     # Returns the position in the vector 'x' where the scalar 'Q' is located
     Qposition <- function(x, Q) {     
       Q.dist  <- abs(x - Q)
       Q.index <- which.min( Q.dist ) 
       return(Q.index)       
     } # end

     # index of those elements that are present both in 'x' and 'y' (NON- NA values)
     vi <- valindex(sim, obs)
     
     # Filtering 'obs' and 'sim', selecting only those pairs of elements 
	 # that are present both in 'x' and 'y' (NON- NA values)
     obs <- as.numeric(obs[vi])
     sim <- as.numeric(sim[vi])
     
     require(hydroTSM) # for using the 'fdc' function
		      
     # Computing the FDC for simulations and observations
     obs.fdc <- fdc(obs, plot=FALSE)
     sim.fdc <- fdc(sim, plot=FALSE)
     
     # Finding the flow value corresponding to the 'lQ.thr' pbb of excedence
     obs.lQ <- obs[Qposition(obs.fdc, lQ.thr)]
     obs.hQ <- obs[Qposition(obs.fdc, hQ.thr)]
     
     sim.lQ <- sim[Qposition(sim.fdc, lQ.thr)]
     sim.hQ <- sim[Qposition(sim.fdc, hQ.thr)]     
     
     denominator <- ( log(obs.hQ) -  log(obs.lQ) )
     
     if ( denominator > 0 ) {
     
       pbiasfdc <- 100 * ( ( ( log(sim.hQ) -  log(sim.lQ) ) / denominator ) - 1 )
     
     } else stop("'log(obs.hQ) -  log(obs.lQ) = 0', it is not possible to compute 'pbiasfdc'") 
     
      if (plot) {
        tmp <- as.matrix(cbind(obs, sim))
        fdc(tmp, lQ.thr=lQ.thr, hQ.thr=hQ.thr, verbose=verbose, ...)
        legend("bottomleft", legend=paste("BiasFDCms=", round(pbiasfdc,1), "%", sep=""), bty="n")
      } # IF end 
     
     return( pbiasfdc )
     
} # 'pbiasfdc.default' end
  
  
pbiasfdc.matrix <- function (sim, obs, lQ.thr=0.7, hQ.thr=0.2, na.rm=TRUE, plot=TRUE, verbose=FALSE, ...){

    pbiasfdc <- rep(NA, ncol(obs))       
          
    pbiasfdc <- sapply(1:ncol(obs), function(i,x,y) { 
                 pbiasfdc[i] <- pbiasfdc.default( x[,i], y[,i], lQ.thr=lQ.thr, hQ.thr=hQ.thr, na.rm=na.rm, plot=plot, verbose=verbose, ...)
            }, x=sim, y=obs )            
           
    return(pbiasfdc)  
     
  } # 'pbiasfdc.matrix' end
  

pbiasfdc.data.frame <- function (sim, obs, lQ.thr=0.7, hQ.thr=0.2, na.rm=TRUE, plot=TRUE, verbose=FALSE, ...){ 
 
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
   
  pbiasfdc.matrix(sim, obs, lQ.thr=lQ.thr, hQ.thr=hQ.thr, na.rm=na.rm, plot=plot, verbose=verbose, ...)
     
} # 'pbiasfdc.data.frame' end 
