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
  