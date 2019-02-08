updateParameter = function(Data){
  
issum = matrix(data = 0, numbin, 1)
edges  = c(seq(0, maxdist, binwidth),Inf)
alphaj = array(0, dim = c(6,numbin))
betaj = array(0, dim = c(6,numbin))

#collect in distance bins
Data$RBin = ceiling(Data$R / binwidth) * binwidth



# summing up the Is for each distance bin
for(j in 1:numbin) {
  issum[j] = sum(Data$Is[Data$RBin == j * binwidth])
}

# number of observation in each bin
n = hist(Data$RBin, breaks = edges, plot = FALSE)$counts[1:numbin]


# posterior estimates
alphaj = alphaj0 + issum
betaj = betaj0 + I0*n - issum


pjHat = (alphaj) / (alphaj + betaj)



# choose only the updated p's ???
pjHat = pjHat[issum != 0]
dfit = distance[issum != 0]

fitData = data.frame(pjHat= pjHat, R = dfit)

return(list(alphaj, betaj, fitData))
}