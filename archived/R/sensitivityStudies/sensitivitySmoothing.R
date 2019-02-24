# SENSITIVITY STUDY
# smoothing function

# Initialization####

source("/functions/initialization.R")

packageList = c("caret", "beepr","nls2","parallel","foreach","doParallel")

initialization(packageList)

# Initializing variables ####
I0 = 9
maxdist = 300
binwidth = 10
numbin = maxdist/binwidth

distance = seq(binwidth/2, maxdist, binwidth)


alphaj0 = rep(1, length(distance))
betaj0 = rep(1, length(distance))

formulaVector = c("pjHat ~ (gamma1/(gamma1+R))^gamma2",
                 "pjHat ~ (gamma1/R)^gamma2",
                 "pjHat ~ gamma1*R + gamma2",
                 "pjHat ~ gamma1*exp(gamma3*R)")

numberData = seq(10, 700, 10)

#number of folds
k = 10 

# number of permutations
perms = 1000

#Progress Bar
#pb = winProgressBar(title="Data Size Binsize", 
#                    label="0% done", 
#                    min=0, 
#                    max=100, 
#                    initial=0)

total = length(numberData)*length(formulaVector)*k*perms

# Computation ####
load("data/syntheticData.Rdata")

no_cores = detectCores() - 1

# Initiate cluster
cl = makeCluster(no_cores)
registerDoParallel(cl)

#start time
strt=Sys.time()

system.time({
  #setup parallel backend to use 8 processors
  
  
  result4=foreach (dataSizeIter = 1:length(numberData),
                   .packages = packageList,
                   .combine = "rbind") %dopar%{
                     
                     result3 = foreach (formulaIter = 1:length(formulaVector),
                                        .packages = packageList,
                                        .combine = "rbind") %dopar%{
                                          
                                          result2 = foreach(permIter = 1:perms,
                                                            .packages = packageList,
                                                            .combine = "rbind") %dopar%{
                                                              
                                                              currentData = Data[sample(1:nrow(Data),nrow(Data)),]
                                                              currentData = currentData[1:numberData[dataSizeIter],]
                                                              
                                                              
                                                              inTraining = createFolds(currentData$R, k = 10, list = T)
                                                              
                                                              result1 = foreach (cvIter = 1:k,
                                                                                 .packages = packageList,
                                                                                 .combine = "rbind") %dopar%{
                                                                                   
                                                                                   #partial = cvIter + (permIter-1)*k + (binSizeIter-1)*k*perms + (dataSizeIter-1)*k*perms*length(binsize)
                                                                                   #info = sprintf("%1.6f%% done", round((partial/total)*100, digits = 6)) 
                                                                                   #setWinProgressBar(pb, (partial/total)*100, label=info) 
                                                                                   
                                                                                   
                                                                                   trainData = currentData[ -inTraining[[cvIter]],]
                                                                                   testData = currentData[inTraining[[cvIter]],]
                                                                                   
                                                                                   parameters = updateParameter(trainData)
                                                                                   #formulae =  pjHat ~ (gamma1/(gamma1+R))^gamma2
                                                                                   
                                                                                   formulae = formulaVector[formulaIter]
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
#Time difference of 6.104916 hours

#close(pb)
beep()

save(meanDiffTrain, meanDiffTest, file = "data/smoothing.RData")

# Figures ####

load("data/smoothing.RData")

#pdf(file = "./thesis/figures/smoothing.pdf",
#    height = 9,
#    width = 16,
#    pointsize = 25)


plot(numberData, meanDiffTest[,1], 
     ylim = c(0.6,1), 
     main = "amount of Data",
     ylab = "mean squared error",
     xlab = "# of Data",
     xaxt = "n",
     xaxs="i", 
     yaxs="i",
     pch = 19)
points(numberData, meanDiffTest[,2], 
       col = "Red",
       pch = 19)
points(numberData, meanDiffTest[,3], 
       col = "Green",
       pch = 19)
points(numberData, meanDiffTest[,4], 
       col = "Blue",
       pch = 19)
lines(numberData, meanDiffTest[,1], col = "Black", lwd = 5)
lines(numberData, meanDiffTest[,2], col = "Red", lwd = 5)
lines(numberData, meanDiffTest[,3], col = "Green", lwd = 5)
lines(numberData, meanDiffTest[,4],col = "Blue", lwd = 5)

axis(1, at=seq(0,max(numberData), by = 10), labels = NA)
axis(lwd=0,side=1,line=-0.4, at=seq(0,max(numberData), by = 20))
#axis(2, at=seq(0,250, by = 20), las = 2)
box()
legend("topright", pch = c(19, 19, 19,19), 
       col = c("red", "black", "blue", "green"), 
       legend = c("gaussian", "uniform", "far", "near"))

#dev.off()


plot(numberData, meanDiffTrain[,1], 
     ylim = c(0.6,1), 
     main = "amount of Data",
     ylab = "mean squared error",
     xlab = "# of Data",
     xaxt = "n",
     xaxs="i", 
     yaxs="i",
     pch = 19)
points(numberData, meanDiffTrain[,2], 
       col = "Red",
       pch = 19)
points(numberData, meanDiffTrain[,3], 
       col = "Green",
       pch = 19)
points(numberData, meanDiffTrain[,4], 
       col = "Blue",
       pch = 19)
lines(numberData, meanDiffTrain[,1], col = "Black", lwd = 5)
lines(numberData, meanDiffTrain[,2], col = "Red", lwd = 5)
lines(numberData, meanDiffTrain[,3], col = "Green", lwd = 5)
lines(numberData, meanDiffTrain[,4],col = "Blue", lwd = 5)

axis(1, at=seq(0,max(numberData), by = 10), labels = NA)
axis(lwd=0,side=1,line=-0.4, at=seq(0,max(numberData), by = 20))
#axis(2, at=seq(0,250, by = 20), las = 2)
box()
legend("topright", pch = c(19, 19, 19,19), 
       col = c("red", "black", "blue", "green"), 
       legend = c("gaussian", "uniform", "far", "near"))
