# SENSITIVITY STUDY
# Prior

# Initialization####

source("functions/initialization.R")

packageList = c("caret", "beepr","nls2","parallel","foreach","doParallel")

initialization(packageList)

# Initializing variables ####
I0 = 9
maxdist = 300


numberData = seq(10, 700, 10)

binwidth = 10
numbin = round(maxdist/binwidth)

distance = seq(binwidth/2, numbin*binwidth, binwidth)

#number of folds
k = 10 

perms = 1000

# Prior

alphaj0Uniform = rep(1, length(distance))
betaj0Uniform = rep(1, length(distance))

# 
alphaj0Tsapanos = (1/(1+distance/3))^(1/I0)
betaj0Tsapanos = 1 - alphaj0Tsapanos
#

priorMean = (I0 -koveslighety(distance,10))/I0
priorVariance = 0.007
alphaj0Koveslighety =  ((1 - priorMean) / priorVariance - 1 / priorMean) * priorMean ^ 2
betaj0Koveslighety = alphaj0Koveslighety * (1 / priorMean - 1)

priorMean = exp(-0.009*distance)
alphaj0Exp =  ((1 - priorMean) / priorVariance - 1 / priorMean) * priorMean ^ 2
betaj0Exp = alphaj0Exp * (1 / priorMean - 1)

priorMean = -1/300*distance + 1
alphaj0Linear =  ((1 - priorMean) / priorVariance - 1 / priorMean) * priorMean ^ 2
betaj0Linear = alphaj0Linear * (1 / priorMean - 1)

alphaArray = cbind(alphaj0Uniform,alphaj0Tsapanos,alphaj0Koveslighety,alphaj0Exp,alphaj0Linear)
betaArray = cbind(betaj0Uniform,betaj0Tsapanos,betaj0Koveslighety,betaj0Exp,betaj0Linear)

#Progress Bar
#pb = winProgressBar(title="Data Size Binsize", 
#                    label="0% done", 
#                    min=0, 
#                    max=100, 
#                    initial=0)

total = length(numberData)*(dim(alphaArray)[2])*k*perms

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
                     
                     result3 = foreach (priorIter = 1:(dim(alphaArray)[2]),
                                        .packages = packageList,
                                        .combine = "rbind") %dopar%{
                                          
                                          
                                          # Prior
                                          alphaj0 = alphaArray[,priorIter]
                                          betaj0 = betaArray[,priorIter]
                                          
                                          
                                          
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
#Time difference of 2.283977 hours

#close(pb)
beep()

save(meanDiffTrain, meanDiffTest, file = "data/prior.RData")


# Figures ####

load("data/prior.RData")

pdf(file = "../thesis/Figures/prior.pdf",
   height = 9,
    width = 16,
    pointsize = 25)

plot(numberData, meanDiffTest[,1], 
     ylim = c(0.75,2), 
     xlim = c(0,700),
     main = "prior",
     ylab = "mean absolute error",
     xlab = "data size",
     #xaxt = "n",
     yaxt = "n",
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
points(numberData, meanDiffTest[,5],
       col = "magenta",
       pch = 19)
lines(numberData, meanDiffTest[,1], col = "Black", lwd = 5)
lines(numberData, meanDiffTest[,2], col = "Red", lwd = 5)
lines(numberData, meanDiffTest[,3], col = "Green", lwd = 5)
lines(numberData, meanDiffTest[,4],col = "Blue", lwd = 5)
lines(numberData, meanDiffTest[,5],col = "Magenta", lwd = 5)

#axis(1, at=seq(0,700, by = 100), labels = seq(0,max(numberData),100))
#axis(1, at=seq(0,max(numberData),50), labels = NA)
axis(2, at=c(0.75,1,1.5,2), las = 2, labels = c(0.75,1,1.5,2))
axis(2, at=seq(0.75,1.75,0.5), labels = NA)

box()
legend("topright", pch = c(19, 19, 19,19, 19), 
       col = c("red", "black", "blue", "green","magenta"), 
       legend = c("Tsapanos", "uniform", "Exponential", "Koveslighety", "Linear"))

dev.off()

# Prior Figures ####



pdf(file = "../thesis/Figures/priorMeans.pdf",
    height = 9,
    width = 16,
    pointsize = 25)

layout(matrix(c(1,2,1,3,1,4), ncol = 3L))
par(mar = c(4,4,4,4),oma = c(0, 0, 2, 0))


matplot(seq(10,300,10), means, 
        type = "l",
        lwd = 5, 
        lty = 1,
        ylab = "p",
        xlab = "distance R [km]",
        main = "mean of prior distributions",
        xaxs="i", yaxs="i",
        xaxt = "n",
        yaxt = "n",
        ylim = c(0,1),
        xlim = c(10,300),
        col = c("black", "red", "green","blue", "magenta"))
axis(1, at=seq(0,300,50), labels = seq(0,300,50))
axis(1, at=seq(0,300,10), labels = NA)
axis(2, at=seq(0,1,0.2), las = 2, labels = seq(0,1,0.2))
axis(2, at=seq(0,1,0.1), labels = NA)

matplot(seq(0,1,0.01),sapply(1:5, function(x) dbeta(seq(0,1,0.01),alphaArray[1,x], betaArray[1,x])), 
        type="l",
        lwd = 5, 
        lty = 1,
        main = sprintf("distance: %d km", 10),
        ylab = "density",
        xlab = "p",
        col = c("black", "red", "green","blue", "magenta")
)

matplot(seq(0,1,0.01),sapply(1:5, function(x) dbeta(seq(0,1,0.01),alphaArray[15,x], betaArray[15,x])),
        type="l",
        lwd = 5, 
        lty = 1,
        main = sprintf("distance: %d km", 150),
        ylab = "density",
        xlab = "p",
        col = c("black", "red", "green","blue", "magenta")
)

matplot(seq(0,1,0.01),sapply(1:5, function(x) dbeta(seq(0,1,0.01),alphaArray[30,x], betaArray[30,x])),
        type="l",
        lwd = 5, 
        lty = 1,
        main = sprintf("distance: %d km", 300),
        ylab = "density",
        xlab = "p",
        col = c("black", "red", "green","blue", "magenta")
)

dev.off()

for (i in seq(1,30,1)){
  
  matplot(seq(0,1,0.005),sapply(1:5, function(x) dbeta(seq(0,1,0.005),alphaArray[i,x], betaArray[i,x])), 
          type="l",
          lwd = 5, 
          lty = 1,
          main = sprintf("distance: %d km", i*10),
          ylab = "density",
          xlab = "p",
          col = c("black", "red", "green","blue", "magenta")
          )
  
  Sys.sleep(1)
}



jpeg(file = "../gif/foo%02d.jpg", type = "cairo",pointsize = 12,quality = 100)
for (i in seq(1,30,1)){
  
  matplot(seq(0,1,0.005),sapply(1:5, function(x) dbeta(seq(0,1,0.005),alphaArray[i,x], betaArray[i,x])), 
          type="l",
          lwd = 5, 
          lty = 1,
          main = sprintf("distance: %d km", i*10),
          ylab = "density",
          xlab = "p",
          col = c("black", "red", "green","blue", "magenta")
  )
  
  #Sys.sleep(1)
}
dev.off()

make.mov = function(){
  unlink("plot.gif")
  system("convert -delay 50.0 ../gif/foo*.jpg ../gif/plot.gif")
}
make.mov()
