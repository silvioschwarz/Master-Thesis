#------------------------------------------------------#
#
#   Intensity Distribution
#   Compute binomial Distribution of Intensities
#
#   Silvio Schwarz
#   (C) 2014-2016 Silvio Schwarz. GNU GPL 3.
#------------------------------------------------------#

initialization = function(packageList){
  # Clear the memory 
  rm(list = ls())
  
  #set working Directory
  path = getwd()
  setwd(path)
  
  # read in additional functions
  pathnames = list.files(pattern="[.]R$", path="functions", full.names=TRUE)
  sapply(pathnames, FUN=source)
  
  #load packages
  loadPackages(packageList) 
  
  #load default parameters
  load("data/R.default.par.RData") 
  par(par.defaults)
  
  
  # setting "stylesheet"
  par(pch = 19, ann = T, cex = 1.5, cex.main = 1.5,font.lab = 2, xaxs="i", yaxs="i")
  
}