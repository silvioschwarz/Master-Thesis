# SENSITIVITY STUDY
# all vs all

# Initialization####
path = getwd()
setwd(path)

# Clear the memory 
rm(list = ls())

# read in additional functions
pathnames = list.files(pattern="[.]R$", path="/functions", full.names=TRUE)
sapply(pathnames, FUN=source)

load("R.default.par.RData"); 
par(par.defaults)

# setting "stylesheet"
par(pch = 19, ann = T, cex = 1.5, xaxs="i", yaxs="i")

# Computation ####
I0 = 9
maxDist = 300


# Initializing variables
numberData = seq(10, 1000, 10)
noise = seq(0,2,0.1)
binsize = seq(1,100,1)
priorVar = seq(0.001,0.01,0.0005)

temp = length(noise) * length(numberData) *length(binsize) * length(priorVar)

#fitPost = vector("list", length(noise))
meanDiff = array(rep(NA, temp), c(length(priorVar),length(binsize),length(noise), length(numberData)))
meanDiffTest = array(rep(NA, temp), c(length(priorVar),length(binsize),length(noise), length(numberData)))
meanDiffTrain = array(rep(NA, temp), c(length(priorVar),length(binsize),length(noise), length(numberData)))


coeffPost = matrix(0, nrow = 6, ncol = 2)

load("syntheticData.Rdata")

#set.seed(5555)

#Data = Data[sample(seq(1, nrow(Data), 1), 1000, replace = F),]
for(v in seq(1,length(priorVar),1)){
  for (u in seq(1,length(binsize),1)){
    
    for (t in seq(1,length(noise),1)){ 
      
      for (i in seq(1,length(numberData),1)){
        
        tryCatch({
          maxdist = binsize[u] * floor(maxDist/binsize[u])
          
          binwidth = binsize[u]
          numbin = maxdist/binwidth
          
          distance = seq(binwidth/2, maxdist, binwidth)
          
          # Tsapanos Prior
          
          priorMean = (1/(1+distance/60))^(1/I0)
          
          maxVariance = priorMean*(1-priorMean)
          
          bound1 = priorMean*(1-priorMean)^2/(2-priorMean)
          bound2 = priorMean^2*(1-priorMean)/(1+priorMean)
          
          
          priorVariance = priorVar[v]
          
          
          alphaj0 =  ((1 - priorMean) / priorVariance - 1 / priorMean) * priorMean ^ 2
          betaj0 = alphaj0 * (1 / priorMean - 1)
          
          #alphaj0 = priorMean
          #betaj0 = 1-alphaj0
          
          
          
          currentData = Data[Data$R<= maxdist,]
          currentData = currentData[1:numberData[i],]
          
          #temp = sapply(probability, intensityDistribution, I0= I0)
          currentData$Is = round(currentData$Is + rnorm(nrow(currentData),0, noise[t]))
          
          
          fitPost = posteriorDistribution(currentData)
          predictData = sapply(predict(fitPost, list(dfit = currentData$R)), intensityDistribution, I0 = 9)
          predictData = apply(predictData, 2, which.max)
          meanDiff[v,u,t,i] = mean(abs(predictData - currentData$Is))
          
          # make folds
          
          
          k = 10 #number of folds
          index = sample(1:nrow(currentData))
          folds = split(index, 1:k)
          
          
          meanTrain = rep(NA, k)
          meanTest = rep(NA, k)
          
          for (j in seq(1,k,1)){
            
            trainData = currentData[-folds[[j]],]
            testData = currentData[folds[[j]],]
            
            fitPostCV = posteriorDistribution(trainData)
            
            
            predictTest = sapply(predict(fitPostCV, list(dfit = testData$R)), intensityDistribution, I0 = 9)
            predictTest = apply(predictTest, 2, which.max)
            predictTrain = sapply(predict(fitPostCV, list(dfit = trainData$R)), intensityDistribution, I0 = 9)
            predictTrain = apply(predictTrain, 2, which.max)
            meanTest[j] = mean(abs(predictTest - testData$Is))
            meanTrain[j] = mean(abs(predictTrain - trainData$Is))
          }
        }, error=function(e){cat("ERROR :",conditionMessage(e), "/n")})
        meanDiffTest[v,u,t,i] = mean(meanTest)
        meanDiffTrain[v,u,t,i] = mean(meanTrain)
      }
    }
  }
}
save(meanDiff, meanDiffTrain, meanDiffTest, file = "allVsAll.RData")



# figures ####
arrayInd(which.min(meanDiff), dim(meanDiff))

load("allVsAll.RData")
loadPackages("fields")
numberData = seq(10,500, 1)


par(mar=c(4,4,4,5))
image(t(meanDiffTest[1,1,,]),
      main = "Bin size vs Data size",
      xlab = "amount of data",
      ylab= "binsize [km]",
      axes=F,
      col='transparent')
axis(1, at=seq(0,1, by = 0.2), labels = seq(0,max(numberData),200))
axis(2, at=seq(0,1, by = 0.25), las = 2, labels = seq(0,max(binsize),25))

image.plot(t(meanDiffTest[1,1,,]),add=T,legend.mar=3.1,legend.args=list(text="MAE"))

box()



numberData = seq(10, 1000, 10)
noise = seq(0,5,0.1)
binsize = seq(1,100,1)
priorVar = seq(0.1,1,0.05)


for(t in seq(1,length(binsize),1)){ 
  
  
  #gap = meanDiffTest-meanDiffTrain
  
  plot(numberData, meanDiffTrain[2,t,1,] , 
       #ylim = c(-0.05,1), 
       main = "amount of Data",
       ylab = "mean squared error",
       xlab = "# of Data",
       xaxt = "n",
       pch = 19)
  points(numberData, meanDiffTest[2,t,1,], 
         col = "Red",
         pch = 19)
  #    points(numberData, gap[t,], 
  #           col = "Blue",
  #           pch = 19)
  
  axis(1, at=seq(0,max(numberData), by = 10), labels = NA)
  axis(lwd=0,side=1,line=-0.4, at=seq(0,max(numberData), by = 20))
  #axis(2, at=seq(0,250, by = 20), las = 2)
  box()
  #legend("topright", pch = c(19, 19), 
  #       col = c("red", "black"), 
  #       legend = c("crossvalidation error", "fit error"))
  Sys.sleep(1)
}

dev.off()

# Tables ####
loadPackages("xtable")
tab = rbind(numberData, meanDiffTrain, meanDiffTest)
row.names(tab) = NULL
xtable(tab)

