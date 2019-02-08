# SENSITIVITY STUDY
# point estimator

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
maxdist = 300
binwidth = 10
numbin = maxdist/binwidth

distance = seq(binwidth/2, maxdist, binwidth)

# Tsapanos Prior

priorMean = (1/(1+distance/3))^(1/I0)
priorVariance = 0.01

alphaj0 =  ((1 - priorMean) / priorVariance - 1 / priorMean) * priorMean ^ 2
betaj0 = alphaj0 * (1 / priorMean - 1)


# Initializing variables
numberData = seq(10, 500, 1)

fitErrorPost =  matrix(0, nrow = 6, ncol = 1)
fitPost = vector("list", length(numberData))
meanDiffTest = rep(NA, length(numberData))
meanDiffTrain = rep(NA, length(numberData))
meanDiffTestRound = rep(NA, length(numberData))
meanDiffTrainRound = rep(NA, length(numberData))
modeDiffTest = rep(NA, length(numberData))
modeDiffTrain = rep(NA, length(numberData))
randomDiffTest = rep(NA, length(numberData))
randomDiffTrain = rep(NA, length(numberData))


coeffPost = matrix(0, nrow = 6, ncol = 2)

load("syntheticData.Rdata")

set.seed(5555)
ind = sample(seq(1, nrow(Data), 1), 1000, replace = F)
Data = Data[ind,]

#Data = sample(Data)
for (i in seq(1,length(numberData),1)){
  
  currentData = Data[1:numberData[i],]
  
  ## make folds
  
  
  k = 10 #number of folds
  index = 1:nrow(currentData)
  folds = split(index, 1:k)
  
  # only validation so fold = 1
  meanTrainMean = rep(NA, k)
  meanTestMean = rep(NA, k)
  meanTrainMeanRound = rep(NA, k)
  meanTestMeanRound = rep(NA, k)
  meanTrainMode = rep(NA, k)
  meanTestMode = rep(NA, k)
  meanTrainRandom = rep(NA, k)
  meanTestRandom = rep(NA, k)
  
  for (j in seq(1,k,1)){
    
    trainData = currentData[-folds[[j]],]
    testData = currentData[folds[[j]],]
    
    fitPostCV = posteriorDistribution(trainData)
    
    predictTest = sapply(predict(fitPostCV, list(dfit = testData$R)), intensityDistribution, I0 = 9)
    predictTestMode = apply(predictTest, 2, which.max)
    predictTestMean = colSums(predictTest * matrix(rep(seq(1,9,1),ncol(predictTest)), c(9,ncol(predictTest))))
    predictTestMeanRound = round(colSums(predictTest * matrix(rep(seq(1,9,1),ncol(predictTest)), c(9,ncol(predictTest)))))
    predictTestRandom = sample(1:9, ncol(predictTest), replace=T)
  
    predictTrain = sapply(predict(fitPostCV, list(dfit = trainData$R)), intensityDistribution, I0 = 9)
    predictTrainMode = apply(predictTrain, 2, which.max)
    predictTrainMean = colSums(predictTrain * matrix(rep(seq(1,9,1), ncol(predictTrain)), c(9,ncol(predictTrain))))
    predictTrainMeanRound = round(colSums(predictTrain * matrix(rep(seq(1,9,1), ncol(predictTrain)), c(9,ncol(predictTrain)))))
    predictTrainRandom = sample(1:9, ncol(predictTrain), replace=T)
    
    meanTestMean[j] = mean((predictTestMean - testData$Is)^2)
    meanTrainMean[j] = mean((predictTrainMean - trainData$Is)^2)
    meanTestMeanRound[j] = mean((predictTestMeanRound - testData$Is)^2)
    meanTrainMeanRound[j] = mean((predictTrainMeanRound - trainData$Is)^2)
    meanTestMode[j] = mean((predictTestMode - testData$Is)^2)
    meanTrainMode[j] = mean((predictTrainMode - trainData$Is)^2)
    meanTestRandom[j] = mean((predictTestRandom - testData$Is)^2)
    meanTrainRandom[j] = mean((predictTrainRandom - trainData$Is)^2)
  }
  meanDiffTest[i] = mean(meanTestMean)
  meanDiffTrain[i] = mean(meanTrainMean)
  meanDiffTestRound[i] = mean(meanTestMeanRound)
  meanDiffTrainRound[i] = mean(meanTrainMeanRound)
  modeDiffTest[i] = mean(meanTestMode)
  modeDiffTrain[i] = mean(meanTrainMode)
  randomDiffTest[i] = mean(meanTestRandom)
  randomDiffTrain[i] = mean(meanTrainRandom)
}

# figures ####

pdf(file = "./thesis/figures/meanModeRandom.pdf",
    height = 9,
    width = 16)

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

gap = meanDiffTest-meanDiffTrain

plot(numberData, meanDiffTrain , 
     ylim = c(0,0.12), 
     main = "Mean",
     ylab = "mean squared error",
     xlab = "# of Data",
     xaxt = "n",
     pch = 19)
points(numberData, meanDiffTest, 
       col = "Red",
       pch = 19)
points(numberData, gap, 
       col = "Blue",
       pch = 19)
axis(1, at=seq(0,max(numberData), by = 10), labels = NA)
axis(lwd=0,side=1,line=-0.4, at=seq(0,max(numberData), by = 20))
#axis(2, at=seq(0,250, by = 20), las = 2)
box()
legend("topright", pch = c(19, 19, 19), 
       col = c("red", "black", "blue"), 
       legend = c("crossvalidation error", "fit error", "difference"))

gapRound = meanDiffTestRound-meanDiffTrainRound

plot(numberData, meanDiffTrainRound , 
     ylim = c(0.0,0.15), 
     main = "Mean rounded",
     ylab = "mean squared error",
     xlab = "# of Data",
     xaxt = "n",
     pch = 19)
points(numberData, meanDiffTestRound, 
       col = "Red",
       pch = 19)
points(numberData, gapRound, 
       col = "Blue",
       pch = 19)
axis(1, at=seq(0,max(numberData), by = 10), labels = NA)
axis(lwd=0,side=1,line=-0.4, at=seq(0,max(numberData), by = 20))
#axis(2, at=seq(0,250, by = 20), las = 2)
box()
legend("topright", pch = c(19, 19, 19), 
       col = c("red", "black", "blue"), 
       legend = c("crossvalidation error", "fit error", "difference"))

modeGap = modeDiffTest-modeDiffTrain

plot(numberData, modeDiffTrain , 
     ylim = c(0.0,0.3), 
     main = "Mode",
     ylab = "mean squared error",
     xlab = "# of Data",
     xaxt = "n",
     pch = 19)
points(numberData, modeDiffTest, 
       col = "Red",
       pch = 19)
points(numberData, modeGap, 
       col = "Blue",
       pch = 19)

axis(1, at=seq(0,max(numberData), by = 10), labels = NA)
axis(lwd=0,side=1,line=-0.4, at=seq(0,max(numberData), by = 20))
#axis(2, at=seq(0,250, by = 20), las = 2)
box()
legend("topright", pch = c(19, 19, 19), 
       col = c("red", "black", "blue"), 
       legend = c("crossvalidation error", "fit error", "difference"))

randomGap = randomDiffTest-randomDiffTrain

plot(numberData, randomDiffTrain , 
     ylim = c(0.0,18), 
     main = "Random",
     ylab = "mean squared error",
     xlab = "# of Data",
     xaxt = "n",
     pch = 19)
points(numberData, randomDiffTest, 
       col = "Red",
       pch = 19)
points(numberData, randomGap, 
       col = "Blue",
       pch = 19)

axis(1, at=seq(0,max(numberData), by = 10), labels = NA)
axis(lwd=0,side=1,line=-0.4, at=seq(0,max(numberData), by = 20))
#axis(2, at=seq(0,250, by = 20), las = 2)
box()
legend("topright", pch = c(19, 19, 19), 
       col = c("red", "black", "blue"), 
       legend = c("crossvalidation error", "fit error", "difference"))

dev.off()

# Tables ####
loadPackages("xtable")
tab = rbind(numberData, meanDiffTrain, meanDiffTest)
row.names(tab) = NULL
xtable(tab)

