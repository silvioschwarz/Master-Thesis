# SENSITIVITY STUDY
# Discretization+NumberData

# Initialization####

source("functions/initialization.R")

packageList = c("caret", "beepr","nls2","parallel","foreach","doParallel")

initialization(packageList)

# Initializing variables ####
I0 = 9
maxdist = 300

binsize = seq(1,100,1)
numberData = seq(10, 700, 10)

#number of folds
k = 10 

perms = 1000

#Progress Bar
#pb = winProgressBar(title="Data Size Binsize", 
#                    label="0% done", 
#                    min=0, 
#                    max=100, 
#                    initial=0)

total = length(numberData)*length(binsize)*k*perms

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
  
result3 = foreach (binSizeIter = 1:length(binsize),
                   .packages = packageList,
                   .combine = "rbind") %dopar%{
  
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
#Time difference of 1.858051 days

#close(pb)
beep()

save(meanDiffTrain, meanDiffTest, file = "data/numberDataBinsize.RData")

# Figures ####
load("data/numberDataBinsize.RData")

loadPackages("fields")

pdf(file = "../thesis/Figures/binsizeDatasize.pdf",
    height = 9,
    width = 16,
    pointsize = 25)

par(mar=c(4,4,4,5))
x = seq(30,700,10)
y = seq(1,100,1)
image(meanDiffTest[3:70,],
      main = "Bin size vs Data size",
      xlab = "data size",
      ylab = "bin size [km]",
      axes = F,
      col='transparent')
axis(1, at=c(0,7/67,17/67,27/67,37/67, 47/67, 57/67,1), las = 1, labels = c(30,100,200,300,400,500,600,700)) 
axis(2, at=seq(0,1, by = 0.2), las = 2, labels = c(1,20,40,60,80,100))

image.plot(meanDiffTest[3:70,],add=T,legend.mar=4.1,legend.args=list(text="MAE"))

box()

dev.off()
