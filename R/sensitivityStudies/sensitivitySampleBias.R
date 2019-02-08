# SENSITIVITY STUDY
# Sample Bias

# Initialization ####

source("functions/initialization.R")

packageList = c("caret", "beepr","nls2","parallel","foreach","doParallel")

initialization(packageList)

# Initializing variables ####
I0 = 9
maxdist = 300
binwidth = 10
numbin = maxdist/binwidth

distance = seq(binwidth/2, maxdist, binwidth)

# Prior

alphaj0 = rep(1, length(distance))
betaj0 = rep(1, length(distance))


numberData = seq(10, 700, 10)
dataArray = c("data/syntheticData.RData",
              "data/syntheticDataNormal.RData", 
              "data/syntheticDataLog1.RData",
              "data/syntheticDataLog2.RData")

meanDiffTest = array(rep(NA, length(numberData) * length(dataArray)), c(length(numberData),length(dataArray)))
meanDiffTrain = array(rep(NA, length(numberData) * length(dataArray)), c(length(numberData),length(dataArray)))

alphaj0 = rep(1, length(distance))
betaj0 = rep(1, length(distance))


k = 10 #number of folds


perms = 1000

total = length(numberData)*k*perms*length(dataArray)

#Progress Bar
#pb = winProgressBar(title="Sample Bias", 
#                    label="0% done", 
#                    min=0, 
#                    max=100, 
#                    initial=0)

# Computation ####

#setup parallel backend to use 8 processors
no_cores = detectCores() - 1

# Initiate cluster
cl = makeCluster(no_cores)
registerDoParallel(cl)

#start time
strt=Sys.time()


system.time({
result4 = foreach (distIter = 1:length(dataArray),
                   .packages = packageList,
                   .combine = "rbind")%dopar%{
  
  load(dataArray[distIter])
  
  result3 = foreach (dataSizeIter = 1:length(numberData),
                     .packages = packageList,
                     .combine = "rbind")%dopar%{
  
    meanTest = array(rep(NA, perms*k), c(perms,k ))
    meanTrain = array(rep(NA, perms*k), c(perms, k))
    
    result2 = foreach(permIter = 1:perms,
                      .packages = packageList,
                      .combine = "rbind")%dopar%{
      currentData = Data[sample(1:nrow(Data),nrow(Data)),]
      currentData = currentData[1:numberData[dataSizeIter],]
      
      
      inTraining = createFolds(currentData$R, k = 10, list = T)
      result1 = foreach (cvIter = 1:k,
                         .packages = packageList,
                         .combine = "rbind")%dopar%{
        
        #partial = cvIter + (permIter-1)*k + (dataSizeIter-1)*k*perms + (distIter-1)*k*perms*length(numberData)
        #info = sprintf("%1.3f%% done", round((partial/total)*100, digits = 3)) 
        #setWinProgressBar(pb, (partial/total)*100, label=info) 
        
        
        trainData = currentData[ -inTraining[[cvIter]],]
        testData = currentData[inTraining[[cvIter]],]
        
        parameters = updateParameter(trainData)
        formulae =  "pjHat ~ (gamma1/(gamma1+R))^gamma2"        
        fitPostCV = posteriorDistribution(parameters[[3]], formulae)
          
          #print(fitPostCV)
          testProbabilities = sapply(predict(fitPostCV, list(R = testData$R)), intensityDistribution, I0 = 9)
          
          #mode
          predictTest = apply(testProbabilities, 2, which.max)
          #mean
          #predictTest = round(colSums(testProbabilities * matrix(rep(seq(1,9,1), ncol(testProbabilities)), c(9,ncol(testProbabilities)))))
          
          
          trainProbabilities = sapply(predict(fitPostCV, list(R = trainData$R)), intensityDistribution, I0 = 9)
          #mode
          predictTrain = apply(trainProbabilities, 2, which.max)
          #mean
          #predictTrain = round(colSums(trainProbabilities * matrix(rep(seq(1,9,1), ncol(trainProbabilities)), c(9,ncol(trainProbabilities)))))
          
          
          meanTest = mean(abs(predictTest-testData$Is), na.rm = T)
          meanTrain = mean(abs(predictTrain-trainData$Is), na.rm =T)
          
          c(meanTest, meanTrain)
                         }
      c(mean(result1[,1], na.rm = T),mean(result1[,2], na.rm = T))
                      }
    c(mean(result2[,1], na.rm = T),mean(result2[,2], na.rm = T))
    
                     }
  c(result3[,1],result3[,2])
                   }
stopCluster(cl)
})

splitIndex = dim(result4)[2]

meanDiffTest = result4[,1:(splitIndex/2)]
meanDiffTrain = result4[,(splitIndex/2+1):splitIndex]
print(Sys.time()-strt)
#Time difference of 3.368021 hours
#close(pb)
beep()

save( meanDiffTrain, meanDiffTest, file = "data/sampleBias.RData")

# Figures ####


pdf(file = "../thesis/Figures/sampleBias.pdf",
    height = 9,
    width = 16,
    pointsize = 25)

layout(matrix(c(1,2,3,4,5,5,5,5), 2, 4, byrow = TRUE))

edges = seq(0,300,20)

load(dataArray[1])
hist(Data$R, breaks = edges, col = "Black", border = "White", xlab = "R", main = "", xlim = c(0,300), ylim = c(0, 200),xaxs="i", yaxs="i")
load(dataArray[2])
hist(Data$R, breaks = edges, col = "Red", xlab = "R", main = "", xlim = c(0,300), ylim = c(0, 200),xaxs="i", yaxs="i")
load(dataArray[3])
hist(Data$R, breaks = edges, col = "Green", xlab = "R", main = "", xlim = c(0,300), ylim = c(0, 200),xaxs="i", yaxs="i")
load(dataArray[4])
hist(Data$R, breaks = edges, col = "Blue", xlab = "R", main = "", xlim = c(0,300), ylim = c(0, 200),xaxs="i", yaxs="i")

par(mar=c(4,4,2,1))
load("data/sampleBias.RData")
plot(numberData, meanDiffTest[1,], 
     ylim = c(0.75,0.9), 
     xlim = c(0, 700),
     ylab = "mean absolute error",
     xlab = "data size",
     xaxt = "n",
     xaxs="i", 
     yaxs="i",
     pch = 19)
points(numberData, meanDiffTest[2,], 
       col = "Red",
       pch = 19)
points(numberData, meanDiffTest[3,], 
       col = "Green",
       pch = 19)
points(numberData, meanDiffTest[4,], 
       col = "Blue",
       pch = 19)
lines(numberData, meanDiffTest[1,], col = "Black", lwd = 5)
lines(numberData, meanDiffTest[2,], col = "Red", lwd = 5)
lines(numberData, meanDiffTest[3,], col = "Green", lwd = 5)
lines(numberData, meanDiffTest[4,],col = "Blue", lwd = 5)

axis(1, at=seq(0,max(numberData), by = 10), labels = NA)
axis(lwd=0,side=1,line=-0, at=seq(0,max(numberData), by = 50))
#axis(2, at=seq(0,250, by = 20), las = 2)
box()
#legend("topright", pch = c(19, 19, 19,19), 
#       col = c("red", "black", "blue", "green"), 
#       legend = c("gaussian", "uniform", "far", "near"))

dev.off()
