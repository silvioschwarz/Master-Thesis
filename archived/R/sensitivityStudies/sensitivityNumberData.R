# SENSITIVITY STUDY
# Number of data

# Initialization####

source("functions/initialization.R")

packageList = c("caret", "beepr","nls2","parallel","foreach","doParallel")

initialization(packageList)


# Initializing variables###
I0 = 9
maxdist = 300
binwidth = 10
numbin = maxdist/binwidth

distance = seq(binwidth/2, maxdist, binwidth)


alphaj0 = rep(1, length(distance))
betaj0 = rep(1, length(distance))

load("data/syntheticData.RData")


numberData = seq(10, nrow(Data), 10)

meanDiffTest = rep(NA, length(numberData))
meanDiffTrain = rep(NA, length(numberData))
k = 10 #number of folds

meanTrain = rep(NA, k)
meanTest = rep(NA, k)
betaSD = rep(NA, k)

perms = 1000

total = length(numberData)*k*perms
meanTest = array(rep(NA, total), c(length(numberData),k, perms))
meanTrain = array(rep(NA, total), c(length(numberData),k, perms))


#Progress Bar
#pb = winProgressBar(title="Data Size", 
#                   label="0% done", 
#                 min=0, 
#                  max=100, 
#                  initial=0)

# Computation ####

#setup parallel backend to use 8 processors
no_cores = detectCores() - 2

# Initiate cluster
cl = makeCluster(no_cores, outfile="")
registerDoParallel(cl)

#start time
strt=Sys.time()
#pb = txtProgressBar(min = 0, max = total, style = 3)

result3 = foreach(dataSizeIter = 1:length(numberData), 
                  .packages = packageList,
                  .combine = "rbind") %dopar%{
                   
                    setTxtProgressBar(pb, dataSizeIter)
                    result2 = foreach(permIter = 1:perms, 
                                      .packages = packageList,
                                      .combine = "rbind") %dopar%{
                                        currentData = Data[sample(1:nrow(Data),nrow(Data)),]
                                        currentData = currentData[1:numberData[dataSizeIter],]
                                        
                                        
                                        inTraining = createFolds(currentData$R, k = 10, list = T)
                                        
                                        
                                        result1 = foreach(cvIter = 1:k, 
                                                          .packages = packageList,
                                                          .combine = "rbind")%dopar%{
                                                            
                                                            #partial = cvIter + (permIter-1)*k + (dataSizeIter-1)*k*perms
                                                           # info = sprintf("%1.3f%% done", round((partial/total)*100, digits = 3)) 
                                                           # setTkProgressBar(pb, (partial/total)*100, label=info) 
                                                            
                                                            trainData = currentData[ -inTraining[[cvIter]],]
                                                            testData = currentData[inTraining[[cvIter]],]
                                                            
                                                            
                                                            parameters = updateParameter(trainData)
                                                            formulae =  "pjHat ~ (gamma1/(gamma1+R))^gamma2"
                                                            
                                                            fitPostCV = posteriorDistribution(parameters[[3]], formulae)
                                                            
                                                            
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
                                                            
                                                            cat(paste("Starting iteration",Sys.time()-strt,"\n")) 
                                                            meanTest = mean(abs(predictTest-testData$Is), na.rm = T)
                                                            meanTrain = mean(abs(predictTrain-trainData$Is), na.rm =T)
                                                            
                                                            c(meanTest, meanTrain)
                                                          }
                                        c(mean(result1[,1], na.rm = T),mean(result1[,2], na.rm = T))
                                      }
                    c(mean(result2[,1], na.rm = T),mean(result2[,2], na.rm = T))
                  }

stopCluster(cl)
print(Sys.time()-strt)
close(pb)
beep()



save(meanDiffTrain, meanDiffTest, file = "data/numberData.RData")

# Figure ####
load("data/numberData.RData")

pdf(file = "../thesis/figures/numberData.pdf",
         height = 9,
          width = 16,
    pointsize = 25)


par(mar=c(4,4,2,1))

plot(numberData, meanDiffTrain , 
     ylim = c(0,0.5), 
     xlim = c(0, max(numberData)),
     main = "Data Size",
     ylab = "mean absolute error",
     xlab = "Data Size",
     xaxt = "n",
     pch = 19,
     xaxs="i", 
     yaxs="i"
)
points(numberData, meanDiffTest, 
       col = "Red",
       pch = 19)
lines(numberData, meanDiffTrain, col = "Black", lwd = 5)
lines(numberData, meanDiffTest, col = "Red", lwd = 5)

axis(1, at=seq(0,max(numberData)+1, by = 50), labels = NA)
axis(lwd=0,side=1,line=-0.4, at=seq(0,max(numberData)+1, by = 100))
#axis(2, at=seq(0,250, by = 20), las = 2)
box()
legend("topright", pch = c(19, 19), 
       col = c("red", "black"), 
       legend = c("cross-validation error", "training error" ))

dev.off()


#Prior ####
plot(seq(0,1,0.01), dbeta(seq(0,1,0.01),alphaj0[30], betaj0[30]), 
     # main = paste(as.character(temp2[i]), "km"),
     ylab = "density",
     xlab = "probability p",
     pch = 19,
     xaxt = "n",
     type = "l",
     lwd = 5)

# Tables ####
loadPackages("xtable")
tab = rbind(numberData, meanDiffTrain, meanDiffTest)
row.names(tab) = NULL
xtable(tab)

