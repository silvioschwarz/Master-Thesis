# SENSITIVITY STUDY
# Number of data

# Initialization####

source("functions/initialization.R")

packageList = c("caret", "beepr","nls2","parallel","foreach","doParallel")

initialization(packageList)

# Initializing variables####
I0 = 9
maxdist = 300
binwidth = 10
numbin = maxdist/binwidth

distance = seq(binwidth/2, maxdist, binwidth)

# Uniform Prior

alphaj0 = rep(1, length(distance))
betaj0 = rep(1, length(distance))


noise = seq(0,20,0.1)

meanDiffTest = rep(NA, length(noise))
meanDiffTrain = rep(NA, length(noise))


k = 10 #number of folds

meanTrain = rep(NA, k)
meanTest = rep(NA, k)
betaSD = rep(NA, k)

perms = 1000

#Progress Bar
#pb = winProgressBar(title="Noise", 
#                    label="0% done", 
#                    min=0, 
#                    max=100, 
#                    initial=0)

total = length(noise)*k*perms
meanTest = array(rep(NA, total), c(length(noise),k, perms))
meanTrain = array(rep(NA, total), c(length(noise),k, perms))

# Computation ####
load("data/syntheticData.Rdata")

#setup parallel backend to use 8 processors
no_cores = detectCores() - 1

# Initiate cluster
cl = makeCluster(no_cores)
registerDoParallel(cl)

#start time
strt=Sys.time()


result3 = foreach(noiseIter = 1:length(noise),
                  .packages = packageList,
                  .combine = "rbind")%dopar%{
  
  result2 = foreach(permIter = 1:perms,
                    .packages = packageList,
                    .combine = "rbind")%dopar%{
    Data = Data[sample(1:nrow(Data),nrow(Data)),]
    currentData = Data
    
    currentData$Is = currentData$Is + round(rnorm(nrow(currentData),0, noise[noiseIter]))
    currentData$Is[currentData$Is > I0] = I0
    currentData$Is[currentData$Is < 1] = 1
    
    
    inTraining = createFolds(currentData$R, k = 10, list = T)
    
    result1 = foreach(cvIter = 1:k,
                     .packages = packageList,
                     .combine = "rbind")%dopar%{
      
      #partial = cvIter + (permIter-1)*k + (noiseIter-1)*k*perms
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
#close(pb)
beep()
#Time difference of 48.663 mins


meanDiffTest = result3[,1]
meanDiffTrain = result3[,2]


save( meanDiffTrain, meanDiffTest, file = "data/noise.RData")

# Figures ####
noise = seq(0,20,0.1) + 1
load("data/noise.RData")

pdf(file = "../thesis/figures/noise.pdf",
    height = 9,
    width = 16,
    pointsize = 25)

plot(noise, meanDiffTest, 
     ylim = c(0, 4),
     xlim = c(1, 21),
     main = "Quality of Data",
     ylab = "mean absolute error",
     xlab = "Standard Deviation of Data",
     xaxt = "n",
     type = "b",
     pch = 19,
     xaxs="i", 
     yaxs="i")
lines(noise, meanDiffTest, lwd = 5)

axis(1, at=seq(0,max(noise)+1, by = 1), labels = NA)
axis(lwd=0,side=1,line=-0.4, at=seq(1,max(noise)+1, by = 1))
#axis(2, at=seq(0,250, by = 20), las = 2)
box()

dev.off()


