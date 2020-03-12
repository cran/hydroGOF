## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("hydroGOF")

## ----message=FALSE------------------------------------------------------------
library(hydroGOF)

## -----------------------------------------------------------------------------
require(zoo)
data(EgaEnEstellaQts)
obs <- EgaEnEstellaQts

## -----------------------------------------------------------------------------
sim <- obs 

## -----------------------------------------------------------------------------
gof(sim=sim, obs=obs)

## -----------------------------------------------------------------------------
sim[1:2000] <- obs[1:2000] + rnorm(2000, mean=10)

## ----fig=TRUE, pdf=TRUE, eps=FALSE, fig.width=12, fig.height=8----------------
ggof(sim=sim, obs=obs, ftype="dm", FUN=mean)

## ----fig=TRUE, pdf=TRUE, eps=FALSE, fig.width=12, fig.height=8----------------
ggof(sim=sim, obs=obs, ftype="dm", FUN=mean, cal.ini="1963-01-01")

## -----------------------------------------------------------------------------
sim <- window(sim, start=as.Date("1963-01-01"))
obs <- window(obs, start=as.Date("1963-01-01"))

gof(sim, obs)

## -----------------------------------------------------------------------------
r <- sim-obs

## -----------------------------------------------------------------------------
library(hydroTSM)
smry(r) 

## ----fig=TRUE, pdf=TRUE, fig.width=12, fig.height=14--------------------------
# daily, monthly and annual plots, boxplots and histograms
hydroplot(r, FUN=mean)

## ----fig=TRUE, eval=TRUE, pdf=TRUE, eps=FALSE, fig.width=12, fig.height=10----
# daily, monthly and annual plots, boxplots and histograms
hydroplot(r, FUN=mean, pfreq="seasonal")

## ----echo=FALSE---------------------------------------------------------------
sessionInfo()$platform
sessionInfo()$R.version$version.string 
paste("hydroGOF", sessionInfo()$otherPkgs$hydroGOF$Version)

