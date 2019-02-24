# SENSITIVITY STUDY
# Discretization

# Initialization####

source("functions/initialization.R")

packageList = c("caret", "beepr","nls2","parallel","foreach","doParallel")

initialization(packageList)

# Initializing variables ####

I0 = 9

maxdist = 300

binsize = seq(1,100,1)
meanDiffTest = rep(NA, length(binsize))
meanDiffTrain = rep(NA, length(binsize))

k = 10 #number of folds

perms = 1000

total = length(binsize)*k*perms

#Progress Bar
#pb = winProgressBar(title="Descretization", 
#                    label="0% done", 
#                    min=0, 
#                    max=100, 
#                    initial=0)

# Computation ####

load("data/syntheticData.RData")


#setup parallel backend to use 8 processors
no_cores = detectCores() - 1

# Initiate cluster
cl = makeCluster(no_cores)
registerDoParallel(cl)

#start time
strt=Sys.time()

result3 = foreach (binSizeIter = 1:length(binsize),
                   .packages = packageList,
                   .combine = "rbind")%dopar%{
                     
                     binwidth = binsize[binSizeIter]
                     numbin = round(maxdist/binwidth)
                     
                     distance = seq(binwidth/2, numbin*binwidth, binwidth)
                     
                     # Prior
                     alphaj0 = rep(1, length(distance))
                     betaj0 = rep(1, length(distance))
                     
                     meanTest = array(rep(NA, perms*k), c(perms,k ))
                     meanTrain = array(rep(NA, perms*k), c(perms, k))
                     
                     result2 = foreach(permIter = 1:perms,
                                       .packages = packageList,
                                       .combine = "rbind")%dopar%{
                                         currentData = Data[sample(1:nrow(Data),nrow(Data)),]
                                         
                                         
                                         inTraining = createFolds(currentData$R, k = 10, list = T)
                                         
                                         
                                                             result1 = foreach(cvIter = 1:k,
                                                                               .packages = packageList,
                                                                               .combine = "rbind")%dopar%{
                                                             #update the progress bar
                                                             #partial = cvIter + (permIter-1)*k + (binSizeIter-1)*k*perms
                                                             #info = sprintf("%1.3f%% done", round((partial/total)*100, digits = 3)) 
                                                             #setWinProgressBar(pb, (partial/total)*100, label=info) 
                                                             
                                                             
                                                             
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
#close(pb)
beep()
#Time difference of 1.037169 hours

meanDiffTest = result3[,1]
meanDiffTrain = result3[,2]

save( meanDiffTrain, meanDiffTest, file = "data/binsize.RData")

# Figures ####

load("data/binsize.RData")

pdf(file = "../thesis/figures/binsize.pdf",
    height = 9,
    width = 16,
    pointsize = 25)

plot(binsize, meanDiffTest, 
     ylim = c(0.76, 0.8),
     xlim = c(0,max(binsize)),
     main = "Binsize",
     ylab = "mean absolute error",
     xlab = "binsize [km]",
     xaxt = "n",
     pch = 19,
     xaxs="i", 
     yaxs="i")
lines(binsize, meanDiffTest, lwd = 5)
# points(binsize, meanDiffTest, 
#        col = "Red",
#        pch = 19)
# points(binsize, meanDiffTrain, 
#        col = "Blue",
#        pch = 19)

x =  300/ seq(1,100,1)
x == round(x)
ticks = which(x == round(x))

#abline(v = 300/ceiling(x), col = 4 )

axis(1, at=seq(0,max(binsize), by = 5), labels = NA)
axis(lwd=0,side=1,line=-0.4, at=seq(0,max(binsize), by = 10))
#axis(1, at = 300/ceiling(x), labels = NA, line = 0, col = 4, tcl=0.25) 
#axis(lwd=0,side=1,line=-2.4, col.axis = 4, labels = ceiling(x), at = 300/ceiling(x))
#axis(2, at=seq(0,250, by = 20), las = 2)
box()

dev.off()
