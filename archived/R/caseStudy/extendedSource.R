# EXTENDED SOURCE####

# Initialization####

source("functions/initialization.R")

packageList = c("caret", "beepr","nls2","parallel","foreach","doParallel", "fields")

initialization(packageList)

# Computation ####

#setup parallel backend to use 8 processors
no_cores = detectCores() - 2

# Initiate cluster
cl = makeCluster(no_cores)
registerDoParallel(cl)

#start time
strt=Sys.time()

dataArray = dataArray = c("data/CaseStudyRound.RData",
                          "data/CaseStudyFloor.RData", 
                          "data/CaseStudyCeil.RData",
                          "data/CaseStudyCertain.RData",
                          "data/CaseStudyInflated.RData",
                          "data/CaseStudyInflatedUncertain.RData",
                          "data/CaseStudyRandom.RData")

result4 = foreach(dataIter = 1:length(dataArray),
                  .packages = packageList,
                  .combine = "rbind")%dopar%{
                    
                    #rm(dataComp)
                    
                    load(dataArray[dataIter])
                    dataComp$R = dataComp$RExt
                    I0Range = seq(min(dataComp$I0), max(dataComp$I0),1)
                    
                    result3 = foreach(I0Iter = I0Range, 
                                      .packages = packageList,
                                      .combine = "rbind") %dopar%{
                                        
                                        
                                        I0 = I0Iter
                                        currentData = dataComp[dataComp$I0 == I0,]
                                        
                                        maxdist = ceiling(max(currentData$R))
                                        binwidth = 2
                                        numbin = ceiling(maxdist/binwidth)
                                        
                                        distance = seq(binwidth/2,numbin*binwidth, binwidth)
                                        
                                        # set prior parameters alpha and beta
                                        alphaj0 = rep(1, length(distance))
                                        betaj0 = rep(1, length(distance))
                                        
                                        
                                        
                                        
                                        parameters = updateParameter(currentData)
                                        formulae =  "pjHat ~ (gamma1/(gamma1+R))^gamma2"        
                                        fitPostCV = posteriorDistribution(parameters[[3]], formulae)
                                        
                                        testProbabilities = sapply(predict(fitPostCV, list(R = currentData$R)), intensityDistribution, I0 = I0)
                                        
                                        #mode
                                        predictTest = apply(testProbabilities, 2, which.max)
                                        #mean
                                        #predictTest = round(colSums(testProbabilities * matrix(rep(seq(1,9,1), ncol(testProbabilities)), c(9,ncol(testProbabilities)))))
                                        
                                        meanRotondi = mean(abs(predictTest-currentData$Is), na.rm = T)
                                        
                                        #ULLAHet al 2015 model
                                        I =  1.007*currentData$mag - 2.004*log10(currentData$depth)+ 3.298 - 
                                          2.692*0.5*log10(currentData$R^2/currentData$depth^2 +1) -
                                          0.000423 * (sqrt(currentData$R^2+currentData$depth^2)-currentData$depth)
                                        
                                        ullahError = mean(abs(round(I) - currentData$Is), na.rm = T)
                                        
                                        c(meanRotondi, ullahError)
                                      }
                    
                    test = array(NA, c(2,11))
                    test2 = rbind(result3[,1], result3[,2])
                    test[1:2,I0Range] = test2
                    #colnames(test)= I0Range
                    test
                    
                  }
stopCluster(cl)
print(Sys.time()-strt)
#close(pb)
beep()
#Time difference of 6.396561 secs 
meanDiffRotondi = t(result4[seq(1,13,2),])
meanDiffUllah = t(result4[seq(2,14,2),])


save(meanDiffRotondi, meanDiffUllah, file = "data/extendedSourceCaseStudy.RData")

# Tables####
load("data/extendedSourceCaseStudy.RData")

library(xtable)

meanDiffRotondi = as.table(meanDiffRotondi)
colnames(meanDiffRotondi) = c("round","floor","ceiling","Certain", "Inflated", "Inflated Uncertain","Random")
rownames(meanDiffRotondi) = seq(1,11,1)
xtable(meanDiffRotondi, digits=rep(6, ncol(meanDiffRotondi)+1))

meanDiffUllah = as.table(meanDiffUllah)
colnames(meanDiffUllah) = c("round","floor","ceiling","Certain", "Inflated","Inflated Uncertain","Random")
rownames(meanDiffUllah) = seq(1,11,1)
xtable(meanDiffUllah, digits=rep(6, ncol(meanDiffUllah)+1))

#normalize results to create ranking 
meanMatrix = matrix(rep(rowMeans(meanDiffRotondi, na.rm = T), ncol(meanDiffRotondi)),nrow(meanDiffRotondi),ncol(meanDiffRotondi))

meanRotondi = meanDiffRotondi - meanMatrix



minMatrix = matrix(rep(apply(meanRotondi,1,min,na.rm = T),ncol(meanRotondi)),nrow(meanRotondi),ncol(meanRotondi))
meanRotondi = meanRotondi - minMatrix
maxMatrix = matrix(rep(apply(meanRotondi,1,max,na.rm = T),ncol(meanRotondi)),nrow(meanRotondi),ncol(meanRotondi))
meanRotondi = meanRotondi/maxMatrix

jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

image.plot(t(meanRotondi[nrow(meanRotondi):1,]))


# show transformed vs original site location
data = read.csv("data/dataCentralAsia.dat", sep = "\t")


transCoord = matrix(NA,nrow(data),2)
for(i in 1:nrow(data)){
  dataTrans = data[i,]
  temp = ellipse2Circle(c(dataTrans$F11,dataTrans$F12),
                        c(dataTrans$F21,dataTrans$F22), 
                        dataTrans$siteLON,dataTrans$siteLAT)
  transCoord[i,1] = temp[[1]]
  transCoord[i,2] = temp[[2]]
  
}


ind = unique(data$ID)
for(i in 1:length(ind)){
  
  dataPlot = data[data$ID == ind[i],]
  dataPlotTrans = transCoord[data$ID == ind[i],]
  
  F1 = c(dataPlot$F11[1],dataPlot$F12[1])
  F2=  c(dataPlot$F21[1],dataPlot$F22[1])
  
  plot(dataPlot$siteLON,dataPlot$siteLAT, xlim = c(70,80), ylim =c(37,55), cex = 2)
  points(dataPlotTrans[,1],dataPlotTrans[,2], col = "cyan")
  points(dataPlot$eventLON,dataPlot$eventLAT, col = "Red")
  #plot(F1[1],F1[2], col = "Blue", xlim = c(min(c(F1[1],F2[1])),max(c(F1[1],F2[1]))), 
  #     ylim =c(min(c(F1[2],F2[2])),max(c(F1[2],F2[2]))))
  #points(F2[1], F2[2], col = "Blue")
  #points(dataPlot$eventLON,dataPlot$eventLAT, col = "Red")
  lines(rbind(F1,F2), lwd = 5)

  readline(prompt = "")
  #Sys.sleep(1)
}
