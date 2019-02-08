#------------------------------------------------------#
#
#   Posterior Distribution
#   Update Beta prior and compute posterior
#   distribution of macroseismic intensities for given
#   distance intervals
#
# needs: numbin, maxdist, binwidth, Data, alphaj0, betaj0
#
#   Silvio Schwarz
#   (C) 2014-2016 Silvio Schwarz. GNU GPL 3.
#------------------------------------------------------#

posteriorDistribution = function(fitData, formulae){
  
  #splitstring = strsplit(formulae, "")[[1]]
  
  #gammaLength = length(unique(splitstring[which(splitstring =="g") +5]))
  
  
  #grd = data.frame(matrix(rep(c(-100,100),2),2,gammaLength))
  #names(grd) = sapply(1:gammaLength, function(x)paste("gamma",x,sep = ""))
  
  grd = data.frame(gamma1 = c(-1,50), gamma2=c(0.1,0.5))
  
  formulae =  as.formula(formulae)
  
  fitPost = NULL
  tryCatch({
    
    
    fitPost = nls(formulae, 
                    data=fitData,
                    start = list(gamma1 = 10, gamma2 = 0.1),
                    control = list(maxiter=200, warnOnly=TRUE))
  }
  , error=function(e){NULL})
    
    if (is.null(fitPost)){
      fitPost = nls2(formulae,
                 data = fitData,
                 start = grd,
                 algorithm = "grid-search",
                 control = list(maxiter=100))
    }
    
  

  return(fitPost)
}
