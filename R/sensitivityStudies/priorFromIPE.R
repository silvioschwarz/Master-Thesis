# Initialization####

source("functions/initialization.R")

packageList = c("caret", "beepr","nls2","parallel","foreach","doParallel")

initialization(packageList)


# ####
I0 = 9

numberData = 10000

r =  runif(numberData, 0, 300)

data = (I0 - koveslighety(r,10)) + rnorm(numberData,0,1)


Data = data.frame(Is = round(data), R = r)

Data$Is[Data$Is > I0] = I0
Data$Is[Data$Is < 1] = 1


# mean and 3 distributions####
alpha = (I0 - koveslighety(r,10)) 
beta = koveslighety(r,10) 

meanIPE = alpha/(alpha+beta)
sdIPE = (alpha*beta)/((alpha+beta+1)*(alpha+beta)^2)


layout(matrix(c(1,2,1,3,1,4), ncol = 3L))
par(mar = c(4,4,4,4),oma = c(0, 0, 2, 0))


plot(r,meanIPE,
     xlim = c(0,300),
     ylim = c(0, 1),
     main = "mean"
     )
#plot(r,sdIPE)

rPlot = seq(10,300,10)
alphaPlot = (I0 - koveslighety(rPlot,10)) 
betaPlot = koveslighety(rPlot,10) 

meanIPEPlot = alphaPlot/(alphaPlot+betaPlot)
sdIPEPLot = (alphaPlot*betaPlot)/((alphaPlot+betaPlot+1)*(alphaPlot+betaPlot)^2)

plot(seq(0,1,0.005),dbeta(seq(0,1,0.005),alphaPlot[1], betaPlot[1]), 
    type="l",
    lwd = 5, 
    lty = 1,
    main = sprintf("distance: %d km", 10),
    ylab = "density",
    xlab = "p",
    col = "blue",
    ylim = c(0,5),
    xlim = c(0,1)
)
plot(seq(0,1,0.005),dbeta(seq(0,1,0.005),alphaPlot[15], betaPlot[15]), 
    type="l",
    lwd = 5, 
    lty = 1,
    ylab = "density",
    main = sprintf("distance: %d km", 150),
    xlab = "p",
    col = "blue",
    ylim = c(0,5),
    xlim = c(0,1)
)
plot(seq(0,1,0.005),dbeta(seq(0,1,0.005),alphaPlot[30], betaPlot[30]), 
    type="l",
    lwd = 5, 
    lty = 1,main = sprintf("distance: %d km", 300),
    ylab = "density",
    xlab = "p",
    col = "blue",
    ylim = c(0,5),
    xlim = c(0,1)
)

# GIF ####
par(mfrow=c(2,1))   
par(mar = c(4,4,3,4),oma = c(0, 2, 0, 0))
for (i in seq(1,30,1)){
  
  plot(rPlot, meanIPEPlot,
       xlim = c(0,300),
       ylim = c(0,1),
       xlab = "R [km]",
       ylab = "p",
       main = sprintf("distance: %d km", i*10))
  lines(rPlot, meanIPEPlot)
  lines(x = rep(i*10,length(seq(0,1,0.1))), y = seq(0,1,0.1), col= "Blue", lwd = 5)
  lines(x = seq(0,i*10,1), y = rep(meanIPEPlot[i],i*10+1), col= "Red", lwd = 5)
  
  plot(seq(0,1,0.005),dbeta(seq(0,1,0.005),alphaPlot[i], betaPlot[i]), 
          type="l",
          lwd = 5, 
          lty = 1,
          ylab = "density",
          xlab = "p",
          col = "blue",
       ylim = c(0,5),
       xlim = c(0,1)
  )
  lines(x = rep(meanIPEPlot[i],length(seq(0,dbeta(meanIPEPlot[i],alphaPlot[i], betaPlot[i]),0.001))), y = seq(0,dbeta(meanIPEPlot[i],alphaPlot[i], betaPlot[i]),0.001), col= "Red", lwd = 5)
  Sys.sleep(1)
}

