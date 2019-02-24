# SENSITIVITY STUDY
# Number of data + noise

# Initialization####

source("/functions/initialization.R")

packageList = c("caret", "beepr","nls2","parallel","foreach","doParallel")

initialization(packageList)

# Initializing variables####
I0 = 9
maxdist = 300
binwidth = 10
numbin = maxdist/binwidth

distance = seq(binwidth/2, maxdist, binwidth)

# Prior
alphaj0 = rep(1, length(distance))
betaj0 = rep(1, length(distance))


numberData = seq(10, 700, 10)
noise = seq(0,5,0.1)



k = 10 #number of folds


perms = 1000

#Progress Bar
#pb = winProgressBar(title="Data Size Noise", 
#                    label="0% done", 
#                    min=0, 
#                    max=100, 
#                    initial=0)

total = length(numberData)*length(noise)*k*perms

meanDiffTest = array(rep(NA, length(numberData)*length(noise)), c(length(numberData),length(noise)))
meanDiffTrain = array(rep(NA, length(numberData)*length(noise)), c(length(numberData),length(noise)))

# Computation ####
load("data/syntheticData.Rdata")

#setup parallel backend to use 8 processors
no_cores = detectCores() - 1

# Initiate cluster
cl = makeCluster(no_cores)
registerDoParallel(cl)

#start time
strt=Sys.time()

system.time({
result4 = foreach (dataSizeIter = 1:length(numberData),
                   .packages = packageList,
                   .combine = "rbind")%dopar%{  
  result3 = foreach (noiseIter = 1:length(noise),
                     .packages = packageList,
                     .combine = "rbind")%dopar%{ 
    
    meanTest = array(rep(NA, perms*k), c(perms,k ))
    meanTrain = array(rep(NA, perms*k), c(perms, k))
    
    result2 = foreach(permIter = 1:perms,
                      .packages = packageList,
                      .combine = "rbind")%dopar%{
      currentData = Data[sample(1:nrow(Data),nrow(Data)),]
      currentData = currentData[1:numberData[dataSizeIter],]
      
      
      
      currentData$Is = round(currentData$Is + rnorm(nrow(currentData),0, noise[noiseIter]))  
      currentData$Is[currentData$Is > I0] = I0  
      currentData$Is[currentData$Is < 1] = 1  
      
      
      inTraining = createFolds(currentData$R, k = 10, list = T) 
      
      result1 = foreach(cvIter = 1:k,
                        .packages = packageList,
                        .combine = "rbind")%dopar%{
        
        #partial = cvIter + (permIter-1)*k + (noiseIter-1)*k*perms + (dataSizeIter-1)*k*perms*length(noise)
        #info = sprintf("%1.6f%% done", round((partial/total)*100, digits = 6)) 
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
  c(result3[,1],result3[,2])
}
stopCluster(cl)
})

splitIndex = dim(result4)[2]

meanDiffTest = result4[,1:(splitIndex/2)]
meanDiffTrain = result4[,(splitIndex/2+1):splitIndex]
print(Sys.time()-strt)
#Time difference of 1.775971 days

#close(pb)
beep()

save(meanDiffTrain, meanDiffTest, file = "data/numberDataNoise.RData")

# figures ####

noise = seq(0,5,0.1) +1
load("data/numberDataNoise.RData")
loadPackages("fields")
numberData = seq(10, 700, 10)


pdf(file = "../thesis/figures/numberDataNoise.pdf",
    height = 9,
    width = 16,
    pointsize = 25)

par(mar=c(4,4,4,5))
image(meanDiffTest,
      main = "Noise vs Data",
      xlab = "data size",
      ylab= "standard deviation",
      axes=F,
      col='transparent')
axis(1, at=seq(0,1, by = 1/7), labels = seq(0,max(numberData),100))
axis(2, at=seq(0,1, by = 0.20), las = 2, labels = seq(0,5,1))

image.plot(meanDiffTest,add=T,legend.mar=3.1,legend.args=list(text="MAE"))

box()

dev.off()
