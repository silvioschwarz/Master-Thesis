# Rotondi 2004 Marche-Umbria Sequence

# initilization ####

# Set working directory

path = getwd()
setwd(path)

# Clear the memory 
rm(list = ls())

# read in additional functions
pathnames = list.files(pattern="[.]R$", path="/functions", full.names=TRUE)
sapply(pathnames, FUN=source)

#packages = c("fields", "xtable")

#loadPackages(packages)
#par.defaults = par(no.readonly=TRUE); 
#save(par.defaults, file="R.default.par.RData") 
load("R.default.par.RData"); 
par(par.defaults)

# setting "stylesheet"
par(pch = 19, ann = T, cex = 1)

# Loading data ####
# reading in data
dataPrior = read.csv("./data/DZ47.dat", sep = "/t")
dataPost = read.csv("./data/DZ47.dat", sep = "/t")


# number of data comparison with rotondi data
#test = apply(table(cut(dataPrior$I0, seq(12,4,-1), right = T), 
#                   cut(dataPrior$Is, seq(0,12,1),right = T)),2,rev)

# binning distance
binwidth = 10

# calculating the maximum and minimum distance in the dataset
mindist = floor(min(dataPrior$R))
maxdist = 250
dataPrior = dataPrior[dataPrior$R <= maxdist,]

# number of bins given maximum distance and bin width (!should be integer!)
numbin = (maxdist - mindist) / binwidth

# set distance equal to middle of bins
d = seq(mindist + binwidth / 2, numbin * binwidth, binwidth)
#d = seq(binwidth, numbin * binwidth, binwidth)
# define edges of the bins
edges = seq(mindist, maxdist, binwidth)

# compute delta I
dataPrior$deltaI = dataPrior$I0 - dataPrior$Is

# Computation of the Prior ####

# initialize variables
numint = 5
pj = array(0, dim = c(numint, numbin))
fitPrior = vector("list", numint)

coeffPrior = array(0, dim = c(numint+1, 2))

absError = matrix(0, nrow = numint, ncol = numbin)
fitErrorPrior =  matrix(0, nrow = numint, ncol = 1)

coeffPaperPrior = matrix(c(7.116, 8.207, 8.317, 7.372, 7.026, 5.903, 
                           0.160, 0.412, 0.305,  0.318, 0.284, 0.206),nrow = 6, ncol = 2)
# set distance equal to middle of bins
d = seq(binwidth, numbin * binwidth, binwidth)

# choose value of I0 to calculate the prior for 
for (i in 5:9){ 
  
  
  # select observation with I0 = i
  #if(i==10){  
  #  ioall = which(dataPrior$I0 %in% c(9.5, 10, 10.5, 11))
  #  data_io = dataPrior[ioall,]
  #  io0 = which(data_io$I0 %in% c(10,11))
  #  ioelse = which(data_io$I0 %in% c(9.5, 10.5))
  #}
  #else{
  ioall = which(dataPrior$I0 %in% c(i, i+0.5))
  data_io = dataPrior[ioall,]
  io0 = which(data_io$I0 == i)
  ioelse = which(data_io$I0 %in% c(i+0.5))
  #}
  
  N = hist(data_io$R, 
           breaks = edges, plot = F)$counts
  
  temp = data_io[io0,]
  temp2 = data_io[ioelse,]
  
  # number of null decay per distance bin I0 = certain, Is weighted
  nio0 = hist(temp[temp$deltaI <=0,]$R,
              breaks = edges, plot = F)$counts +
    0.5 * hist(temp[temp$deltaI == 0.5,]$R, 
               breaks = edges, plot = F)$counts
  
  # number of null decay per distance bin I0 = uncertain, Is uncertain
  nioelse = hist(temp2[temp2$deltaI <= -0.5,]$R, 
                 breaks = edges, plot = F)$counts +
    0.75 * hist(temp2[temp2$deltaI == 0,]$R, 
                breaks = edges, plot = F)$counts + 
    0.5 * hist(temp2[temp2$deltaI == 0.5,]$R, 
               breaks = edges, plot = F)$counts + 
    0.25 * hist(temp2[temp2$deltaI == 1,]$R, 
                breaks = edges, plot = F)$counts
  
  
  NIO =  nio0 + nioelse
  
  
  pj[i-4,] = (NIO/N)^(1/i)
  pj[pj == 0] = NA
  
  
  p = pj[i-4,]
  dfit= d[which(is.finite(p))] # only distance where there is data
  p = p[is.finite(p)] # only data that != NA
  
  par(mfrow=c(1,1))                
  fitPrior[[i-4]] = nls(p  ~ (c1/dfit)^c2, 
                        data = data.frame(cbind(dfit, p)), 
                        start = list(c1 = 7,c2 = 0.3))
  
  
  absError[i-4,1:length(p)] = abs(predict(fitPrior[[i-4]]) - p)
  fitErrorPrior[i-4] = sqrt(mean((predict(fitPrior[[i-4]]) - p)^2, na.rm = TRUE))
  
  coeffPrior[i-4,] = coef(fitPrior[[i-4]])
  
  # PLOTTING
  
  plot(d,pj[i-4,], 
       ylim=c(0,1), xlim = c(0, 240), 
       col = 'blue',
       xlab = 'epicentral R [km]', ylab = 'p0',
       xaxt = "n")
  predictDistance = seq(1,240,5)
  predictPrior = predict(fitPrior[[i-4]], list(dfit = predictDistance))
  predictPrior[predictPrior> 0.98] = 0.98
  
  lines(predictDistance, predictPrior,
        col = 'blue', 
        xaxt = "n")
  
  predictPaper = (coeffPaperPrior[i-4, 1]/predictDistance)^coeffPaperPrior[i-4,2]
  predictPaper[predictPaper > 0.98] = 0.98
  lines(predictDistance, predictPaper,
        col = 'red', 
        xaxt = "n")
  
  
  axis(1, at=seq(0,250, by = 10), labels = NA)
  axis(lwd=0,side=1,line=-0.4, at=seq(0,250, by = 20))
  #axis(2, at=seq(0,250, by = 20), las = 2)
  box()
  #if(i==10)  
  #  title('epicentral Intensity: 10 + 11')
  #else
  title(paste('epicentral Intensity: ',as.character(i)))
  legend("topright", 
         c('computed','a priori estimates (zero decay)', 
           "paper"),
         lty = c( NA, "solid","solid"), pch = c("o", "", ""),
         col = c("blue", "blue", "red")) 
  
}

maxAbsError = apply(absError, 1, function(x) max(x[x != 0]))
minAbsError = apply(absError, 1, function(x) min(x[x != 0]))
coeffPrior - coeffPaperPrior

testDistance = seq(0.1,250,0.1)

for (i in 5:9){ 
  predictPrior = predict(fitPrior[[i-4]], list(dfit = testDistance))
  predictPaperPrior = (coeffPaperPrior[i-4,1]/testDistance)^coeffPaperPrior[i-4,2]
  test = sum(abs(predictPrior - predictPaperPrior), na.rm = T)
  print(test)
}

# Computation of the Posterior ####
errorPaper = c(0.005992,0.001112,0.001829,0.005237,0.001086,0.000432)

coeffPaperPost = matrix(c(7.116, 8.833, 8.924, 10.258, 9.052, 11.600,
                          0.160, 0.324,  0.298,  0.283,  0.318, 0.251),nrow = 6, ncol = 2)


par(par.defaults)
alphaj0 = matrix(0, nrow = 6, ncol = numbin)
betaj0 = matrix(0, nrow = 6, ncol = numbin)

coeffPaperPrior = matrix(c(7.116, 8.207, 8.317, 7.372, 7.026, 5.903, 
                           0.160, 0.412, 0.305,  0.318, 0.284, 0.206), nrow = 6, ncol = 2)

# inverting mean and varinace to obtain alpha and beta TODO: how to choose variance
for (i in 5:10){
  d = seq(binwidth, numbin*binwidth, binwidth)
  priorMean = (coeffPaperPrior[i-4,1]/d)^coeffPaperPrior[i-4,2]
  priorMean[priorMean > 0.98] = 0.98
  
  bound1 = ((1-priorMean)^2/(2-priorMean))*priorMean
  bound2 = ((1-priorMean)/(1+priorMean))*priorMean^2
  
  maxConstraint = bound2
  maxConstraint[which(bound1 > bound2)] = bound1[which(bound1 > bound2)]
  minConstraint = bound2
  minConstraint[which(bound1 < bound2)] = bound1[which(bound1 < bound2)]
  
  #priorVariance = minConstraint + (maxConstraint - minConstraint)*0.01
  expVar = c(0.1, 0.023, 0.03, 0.038, 0.0075,0.01)
  #priorVariance = 0.08
  #priorVariance = rep(NA, length(priorMean))
  #priorVariance[1] = errorPaper[i-4]
  #priorVariance[2:numbin] = seq(minAbsError, maxAbsError[i-4],length.out = numbin -1)
  priorVariance = seq(errorPaper[i-4], expVar[i-4],length.out = numbin)
  
  
  
  alphaj0[i-4,] =  ((1 - priorMean) / priorVariance - 1 / priorMean) * priorMean ^ 2
  betaj0[i-4,] = alphaj0[i-4,] * (1 / priorMean - 1)
  
  
  #alphaj0[i-4,] = priorMean * ((1 - priorMean) * priorMean - priorVariance) / priorVariance
  #betaj0[i-4,] = (1 - priorMean) * ((1 - priorMean) * priorMean - priorVariance) / priorVariance
  
}


dataPost = read.csv("./data/D47.dat", sep = "/t")

# calculating the maximum and minimum distance in the dataset
dataPost = dataPost[dataPost$R <= 250,]

# initializing matrixes


fitErrorPost =  matrix(0, nrow = 6, ncol = 1)
fitPost = vector("list", 5)
alphaj = array(0, dim = c(6,numbin))
betaj = array(0, dim = c(6,numbin))

coeffPost = matrix(0, nrow = 6, ncol = 2)

issum = matrix(data = 0, numbin, 1)

# 
for (i in 6:9){
  edges  = seq(0, maxdist, binwidth)
  
  ioall_post = which(dataPost$I0 %in% c(i, i+0.5))
  
  data_iopost = dataPost[ioall_post,]
  
  #collect in distance bins
  data_iopost$R = ceiling(data_iopost$R / binwidth) * binwidth
  
  
  x = seq(binwidth,numbin*binwidth,binwidth)
  
  # summing up the Is for each distance bin
  for(j in 1:numbin) {
    issum[j] = sum(data_iopost$Is[data_iopost$R == j * binwidth])
  }
  
  # number of observation in each bin
  n = hist(data_iopost$R, breaks = edges, plot = FALSE)$counts[1:numbin]
  
  # posterior estimates
  alphaj[i-4,] = alphaj0[i-4,] + issum
  betaj[i-4,] = betaj0[i-4,] + i*n - issum
  
  pjHat = (alphaj0[i-4,] + issum) / (alphaj0[i-4,] + betaj0[i-4, ] + i * n) 
  
  
  
  # choose only the updated p's ???
  
  pjHat = pjHat[issum != 0,]
  dfit = x[issum != 0]
  
  pjHat[pjHat > 0.98] = 0.98
  fitPost[[i-4]] = nls(pjHat ~ (gamma1/dfit)^gamma2, 
                       data=data.frame(cbind(dfit, pjHat)),
                       start = list(gamma1 = 7.0, gamma2 = 0.3))
  
  coeffPost[i-4,] = coef(fitPost[[i-4]])
  
  plot(dfit, pjHat, 
       xlim = c(0, 240), 
       ylim = c(0, 1), 
       col = 'black', 
       pch = 0, 
       axes = FALSE,
       xlab = 'epicentral distance [km]', 
       ylab = 'p')
  
  lines(x, predict(fitPost[[i-4]], list(dfit = x)), col = "black")
  lines(x, (coeffPaperPost[i-4,1]/x)^coeffPaperPost[i-4,2], col = "blue", lty = "dashed")
  #points(dfit,p, pch = 19, col = "black")
  #lines(x, (coeffPaperPrior[i-4,1]/x)^coeffPaperPrior[i-4,2], col = "black", lty = "dashed")
  axis(1, at=seq(0,240,length.out = 25), labels = FALSE)
  axis(1, at=seq(0,240,length.out = 13), labels = seq(0,240,20))
  axis(2, at=seq(0,1,length.out = 6), labels = seq(0,1,0.2))
  box()
  title(paste('epicentral Intensity: ',as.character(i)))
  
  legend("topright", 
         c('computed','a priori estimates (zero decay)', 
           "paper"),
         lty = c( NA, "solid","solid"), pch = c("o", "", ""),
         col = c("blue", "blue", "red")) 
  
  
  fitErrorPost[i-4] = sqrt(mean((predict(fitPost[[i-4]]) - pjHat)^2,na.rm = TRUE))
}


coeffPost[2:6,] - coeffPaperPost[2:6,]


testDistance = seq(0,250,0.1)

for (i in 6:9){ 
  predictPost = predict(fitPost[[i-4]], list(dfit = testDistance))
  predictPaperPost = (coeffPaperPost[i-4,1]/testDistance)^coeffPaperPost[i-4,2]
  test = sum(abs(predictPost - predictPaperPost), na.rm = T)
  print(test)
}

# Comparison with real data####

# make them integer
dataPost$Is = floor(dataPost$Is)
dataPost$I0 = floor(dataPost$I0)

# model

binomDistro = function(R,I0){
  # compute probability distribution
  probability = (coeffPost[I0-4,1]/R)^coeffPost[I0-4,2]
  
  probability[probability > 0.98] = 0.98
  
  prob = dbinom(x = seq(1, I0, 1), size = I0, prob = probability)
  mode = which.max(prob)
  meann = round(sum(prob*seq(1, I0, 1)))
  
  return(append(prob,mode))
  
}
coeffPaperPost = matrix(c(7.116, 8.833, 8.924, 10.258, 9.052, 11.600,
                          0.160, 0.324,  0.298,  0.283,  0.318, 0.251),nrow = 6, ncol = 2)


discr = array(NA, dim = c(5))
for(i in 6:10){
  data_val = dataPost[which(dataPost$I0 %in% c(i)),]
  test = sapply(data_val$R, binomDistro, I0 = i)
  hist(floor(data_val$Is) - test[dim(test)[1],])
  discr[i-5] = mean(floor(data_val$Is) - test[dim(test)[1],])
}





