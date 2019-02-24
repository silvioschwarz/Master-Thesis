# EXTENDED SOURCE

# Initialization####

source("functions/initialization.R")

packageList = c("caret", "beepr","nls2","parallel","foreach","doParallel")

initialization(packageList)

# Computation####
I0 = 9
y = runif(10000, -300, 300)
x = runif(10000, -300, 300)


#plot(x,y)

F1 = c(90,0)
F2 = c(-90,50)

epicentre = (F1+F2)/2

RCircle = sqrt((epicentre[1]-x)^2 + (epicentre[2]-y)^2)

r1 = sqrt((F1[1]-x)^2 + (F1[2]-y)^2)
r2 = sqrt((F2[1]-x)^2 + (F2[2]-y)^2)
c = sqrt(sum((F1-F2)^2))/2
a = (r1+r2) /2
b = sqrt(a^2-c^2)

REllipse = (r1+r2) /2 -c#-a


intensity = round(I0 - koveslighety(REllipse,10))



levels = seq(1,9,length.out = 9)
col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                          "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(9)  
colz = col[intensity]  
#   

layout(matrix(c(2, 1), ncol = 2L), widths = c(5, 1))
par(mar = c(4,1,4,4))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0,9), xaxs = "i", yaxs = "i")

rect(0, seq(0,8,1), 1, seq(1,9,1),col=col,border=col) 
mtext("Intensity", side = 3, line = 1)
box()
axis(4, at=seq(0,9, by = 1), labels = NA)
axis(lwd=0,side=4, at=seq(0.5,8.5, by = 1), labels = seq(1,9,1), las = 1)
par(mar = c(4,4,4,1))
plot(x,y, col = colz)
lines(rbind(F1,F2), lwd = 5)

toDegrees = 180/pi

azimuth = asin(abs((F1-F2))[1]/(2*c))*toDegrees 
theta = 90 - azimuth 

x1 = x
y1 = y
x2 = cos(-theta)*x1-sin(-theta)*y1 
y2 = sin(-theta)*x1+cos(-theta)*y1 

x3 = x2 *b/a
y3 = y2

gamma = atan(y3/x3)-atan(y2/x2)+abs(theta)
x4 = cos(gamma)*x3-sin(gamma)*y3
y4 = sin(gamma)*x3+cos(gamma)*y3



RCircleTransformed = sqrt((epicentre[1]-x4)^2 + (epicentre[2]-y4)^2)

Data= data.frame(Is =intensity, REllipse = REllipse, RCircle = RCircleTransformed)

save(Data, file="data/extendedSource.RData")

# Figure  ####
load("data/extendedSource.RData")

pdf(file = "../thesis/Figures/extendedSource.pdf",
    height = 10,
    width = 16,
    pointsize = 25)

layout(matrix(c(2,4,3,5, 1,1), ncol = 3L), widths = c(5, 5,1.2))
par(mar = c(4,1,4,4),oma = c(0, 0, 2, 0))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0,9), xaxs = "i", yaxs = "i")

rect(0, seq(0,8,1), 1, seq(1,9,1),col=col,border=col) 
mtext("Intensity", side = 3, line = 1, cex = 1.5)
box()
axis(4, at=seq(0,9, by = 1), labels = NA)
axis(lwd=0,side=4, at=seq(0.5,8.5, by = 1), labels = seq(1,9,1), las = 1,font = 2, cex.axis = 1.5)
xplot = cbind(x1,x2,x3,x4)
yplot = cbind(y1,y2,y3,y4)

mtext("Transformation Extended Source", outer = TRUE, cex = 1.5)


for (i in 1:4){
par(mar = c(4,4,4,1))
plot(xplot[,i],yplot[,i], 
     col = colz,
     #xlim = c(-300,300),
     ylim = c(-300,300),
     main = letters[i],
     xlab = "x",
     ylab= "y", 
     asp = 1,
     pch = 19,
     font.lab = 2)
lines(rbind(F1,F2), lwd = 5)
}

dev.off()

# Computing the differences ####


# Initialization####

source("functions/initialization.R")

packageList = c("caret", "beepr","nls2","parallel","foreach","doParallel")

initialization(packageList)

# Computation ####
I0 = 9
maxdist = 300
binwidth = 10
numbin = maxdist/binwidth

distance = seq(binwidth/2, maxdist, binwidth)

# Prior
alphaj0 = rep(1, length(distance))
betaj0 = rep(1, length(distance))



load("data/extendedSource.RData")


k = 10 #number of folds


perms = 1000

#Progress Bar
#pb = winProgressBar(title="Data Size Noise", 
#                    label="0% done", 
#                    min=0, 
#                    max=100, 
#                    initial=0)

total = k*perms*2

meanTest = array(rep(NA, perms*k), c(perms,k ))
meanTrain = array(rep(NA, perms*k), c(perms, k))


meanDiffTest = array(rep(NA, 2), 2)
meanDiffTrain = array(rep(NA, 2), 2)

#setup parallel backend to use 8 processors
no_cores = detectCores() - 1

# Initiate cluster
cl = makeCluster(no_cores)
registerDoParallel(cl)

#start time
strt=Sys.time()

result3 = foreach(i = 2:3,
                  .packages = packageList,
                  .combine = "rbind") %dopar%{

    result2 = foreach(permIter =  1:perms,
                      .packages = packageList,
                      .combine = "rbind") %dopar%{
      Data$R = Data[,i]
      currentData = Data[sample(1:nrow(Data),nrow(Data)),]
      
      
      #plot(currentData$R, currentData$Is)
      
      
      inTraining = createFolds(currentData$R, k = 10, list = T) 
      
      result1 = foreach(cvIter = 1:k,
                        .packages = packageList,
                        .combine = "rbind") %dopar%{
        
        #partial = cvIter + (permIter-1)*k + (i-2)*k*perms# + (dataSizeIter-1)*k*perms*length(noise)
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
        
        cat(paste("Starting iteration",Sys.time()-strt,"/n"))
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

#Time difference of 33.70072 mins
meanDiffTest = result3[,1]
meanDiffTrain = result3[,2]


save(meanDiffTrain, meanDiffTest, file = "data/extendedSourceResults.RData")

#meanDiffTrain
#result.1  result.2 
#0.1551224 0.2623517 
# ####
load("data/extendedSourceResults.RData")
meanDiffTrain

