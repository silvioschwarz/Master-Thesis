# SENSITIVITY STUDY
# Number of data

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
par(pch = 19, ann = T, cex = 1, xaxs="i", yaxs="i")

# Data Generation  Uniform Distribution in R ####
#I0 = 9

#numberData = 1000000

#set.seed(12345)
#r =  runif(numberData, 0, 300)

#data = (I0 - koveslighety(r,10))


#probability = data/I0 
#plot(r, probability)


#Data = data.frame(Is = sapply(probability, binomialDistribution, I0 = 9), R = r, I0 = I0)

#plot(Data$R, Data$Is,
#     xlim = c(0,300),
#     ylim = c(1,9),
#     xaxs="i",
#     yaxs="i") 

#save(Data, file = "syntheticData.RData")


# Computation ####
I0 = 9
maxdist = 300
binwidth = 10
numbin = maxdist/binwidth

distance = seq(binwidth/2, maxdist, binwidth)

# Tsapanos Prior

alphaj0 = (1/(1+distance/3))^(1/I0)
betaj0 = 1 - alphaj0

# Initializing variables
numberData = seq(10, 400, 1)

fitErrorPost =  matrix(0, nrow = 6, ncol = 1)
fitPost = vector("list", length(numberData))
meanDiffTest = rep(NA, length(numberData))
meanDiffTrain = rep(NA, length(numberData))


coeffPost = matrix(0, nrow = 6, ncol = 2)

load("syntheticData.Rdata")

Data = Data[sample(seq(1, nrow(Data), 1), 400, replace = F),]

for (i in seq(1,length(numberData),1)){

currentData = Data[1:numberData[i],]

## make folds
set.seed(5555)

k = 10 #number of folds
index = 1:nrow(currentData)
folds = split(index, 1:k)

meanTrain = rep(NA, 10)
meanTest = rep(NA, 10)

for (j in seq(1,k,1)){
  
  trainData = currentData[-folds[[j]],]
  testData = currentData[folds[[j]],]
  
  fitPostCV = posteriorDistribution(trainData)
 
  predictTest = sapply(predict(fitPostCV, list(dfit = testData$R)), binomialDistribution, I0 = 9)
  predictTrain = sapply(predict(fitPostCV, list(dfit = trainData$R)), binomialDistribution, I0 = 9)
  meanTest[j] = mean((predictTest - testData$Is)^2)
  meanTrain[j] = mean((predictTrain - trainData$Is)^2)
}
meanDiffTest[i] = mean(meanTest)
meanDiffTrain[i] = mean(meanTrain)
}

plot(numberData, meanDiffTrain , ylim = c(0, 1))
points(numberData, meanDiffTest, col = "Red")
points(numberData, meanDiffTest-meanDiffTrain, col = "Blue")


# Folds ####

# make folds
set.seed(5555)

k = 10 #number of folds
index = sample(1:nrow(currentData))
folds = split(index, 1:k)

for (i in length(folds)){
fitPostCV = posteriorDistribution(currentData[-folds[[i]],])
}

coeffPost = coef(fitPost)

plot(dfit, pjHat, 
     xlim = c(0, 240), 
     ylim = c(0, 1), 
     col = 'black', 
     pch = 0, 
     axes = FALSE,
     xlab = 'epicentral distance [km]', 
     ylab = 'p')



predictDistance = seq(1,300,1)
predictPost = predict(fitPost, list(dfit = predictDistance))
predictPost[predictPost> 0.98] = 0.98

lines(predictDistance, predictPost,
      col = 'black', 
      xaxt = "n")

predictTarget = (I0-koveslighety(predictDistance, 10))/I0
predictTarget[predictTarget> 0.98] = 0.98

lines(predictDistance, predictTarget, col = "blue", lty = "dashed")
axis(1, at=seq(0,240,length.out = 25), labels = FALSE)
axis(1, at=seq(0,240,length.out = 13), labels = seq(0,240,20))
axis(2, at=seq(0,1,length.out = 6), labels = seq(0,1,0.2))
box()
title(paste('epicentral Intensity: ',as.character(I0)))

legend("topright", 
       c('computed','a priori estimates (zero decay)', 
         "paper"),
       lty = c( NA, "solid","solid"), pch = c("o", "", ""),
       col = c("blue", "blue", "red")) 

# Validation####
#mean(predictTarget-predictPost)
#sum(predictTarget-predictPost)

predictIntensities = sapply(predict(fitPost, list(dfit = r)), binomialDistribution, I0 = 9)

plot(r, predictIntensities, col=rgb(0, 0, 1, 0.5))
points(Data$R, Data$Is,
       col=rgb(1, 0, 0, 0.1))

mean(predictIntensities - Data$Is)

hist(predictIntensities - Data$Is, plot = F, breaks = seq(-1.5,1.5,1))




folds[[1]]
for (i in length(folds)){
  
  predictIntensities = sapply(predict(fitPost, list(dfit = r[-folds[[i]]])), binomialDistribution, I0 = 9)
  
  
  
}

# direct comparison####

fitDirect = nls(probability ~ (gamma1/(gamma1+r))^gamma2, 
              data=data.frame(cbind(r, probability)),
              start = list(gamma1 = 7.0, gamma2 = 0.3))

coeffPost = coef(fitDirect)

plot(r, probability, 
     xlim = c(0, 240), 
     ylim = c(0, 1), 
     col = 'black', 
     pch = 0, 
     axes = FALSE,
     xlab = 'epicentral distance [km]', 
     ylab = 'p')

predictDistance = seq(1,300,1)
predictDirect = predict(fitDirect, list(r = predictDistance))
predictDirect[predictDirect> 0.98] = 0.98

lines(predictDistance, predictDirect,
      col = 'black', 
      xaxt = "n")

predictTarget = (I0-koveslighety(predictDistance, 10))/I0
predictTarget[predictTarget> 0.98] = 0.98

lines(predictDistance, predictTarget, col = "blue", lty = "dashed")
axis(1, at=seq(0,240,length.out = 25), labels = FALSE)
axis(1, at=seq(0,240,length.out = 13), labels = seq(0,240,20))
axis(2, at=seq(0,1,length.out = 6), labels = seq(0,1,0.2))
box()
title(paste('epicentral Intensity: ',as.character(I0)))

legend("topright", 
       c('computed','a priori estimates (zero decay)', 
         "paper"),
       lty = c( NA, "solid","solid"), pch = c("o", "", ""),
       col = c("blue", "blue", "red")) 


predictIntensities = sapply(predict(fitDirect, list(r = r)), binomialDistribution, I0 = 9)

plot(r, predictIntensities, col=rgb(0, 0, 1, 0.5))
points(Data$R, Data$Is,
       col=rgb(1, 0, 0, 0.1))

mean(predictIntensities - Data$Is)

hist(predictIntensities - Data$Is, plot = F, breaks = seq(-1.5,1.5,1))
