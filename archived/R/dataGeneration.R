# Initialization####

source("functions/initialization.R")

packageList = c("caret", "beepr","nls2","parallel","foreach","doParallel")

initialization(packageList)

# Data Generation  Uniform Distribution ####
I0 = 9

numberData = 700

r =  runif(numberData, 0, 300)

data = (I0 - koveslighety(r,10)) + rnorm(numberData,0,1)

Data = data.frame(Is = round(data), R = r)

Data$Is[Data$Is > I0] = I0
Data$Is[Data$Is < 1] = 1

save(Data, file = "data/syntheticData.RData")

# Figure Data####
pdf(file = "../thesis/Figures/dataDistro.pdf",
    height = 9,
    width = 16,
    pointsize = 25)

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mar=c(4,4,2,1))


plot(Data$R, round(Data$Is),
     xlim = c(0,300),
     ylim = c(1,10),     
     main = "synthetic Intensity Data",
     ylab = expression("Intensity at Site I"[" s"]),
     xlab = "Distance R [km]",
     xaxt = "n",
     xaxs="i", 
     yaxs="i",
     pch = 19,
     font.lab = 2) 

axis(1, at=seq(0,max(Data$R), by = 10), labels = NA)
axis(lwd=0,side=1,line=-0.4, at=seq(0,max(Data$R)+1, by = 20))
box()


hist(Data$R, 
     xlim = c(0,300),
     ylim = c(0,100),     
     main = "Distance Distribution",
     ylab = "Frequency",
     xlab = "Distance R [km]",xaxs="i", yaxs="i")
box()

hist(Data$Is, 
     breaks = seq(0.5,9.5,1),
     ylim = c(0,200),     
     main = "Intensity Distribution",
     xlab = "Intensity",
    xaxs="i", yaxs="i",
    xaxt = "n")
axis(lwd=0,side=1,line=-0.4, at=seq(1,I0,1))
box()

dev.off()

# Example####

pdf(file = "../thesis/Figures/example.pdf",
    height = 9,
    width = 10,
    pointsize = 25)

layout(matrix(c(1,1,1,2,2,2,3,4,5), 3, 3, byrow = TRUE))
par(mar = c(2,3,2,2))
par(mgp = c(2,1,0))

exData = round(data + rnorm(length(data),0,0.5))
ind = sample(seq(1,length(exData),1),300, replace = F)
plot(r[ind], exData[ind],
     xlim = c(0,300),
     ylim = c(0,10),     
     main = "synthetic Intensity Data",
     ylab = expression("Intensity at Site I"[" s"]),
     xlab = "",
     xaxt = "n",
     pch = 19,
     lwd = 5,
     xaxs="i", yaxs="i")

abline(v = seq(10,300,10), col = 4)
axis(1, at=seq(0,max(Data$R)+1, by = 10), labels = NA)
axis(lwd=0,side=1,line=-0.4, at=seq(0,max(Data$R)+1, by = 20))
#axis(2, at=seq(0,250, by = 20), las = 2)
box()

distance = seq(5,295,10)

data = (I0 - koveslighety(distance,10))

probability = data/I0 

par(mar = c(4,3,0,2))
plot(distance, probability,
     xlim = c(0,300),
     ylim = c(0,1),     
     #main = "synthetic Intensity Data",
     ylab = "p",
     xlab = "Distance R [km]",
     xaxt = "n",
     pch = 19,
     lwd = 5,
     xaxs="i", yaxs="i")
lines(distance, probability,lwd=5, col = 4)

abline(v = seq(5,295,10), col = 4)
axis(1, at=seq(5,295, by = 10), labels = NA)
axis(lwd=0,side=1,line=-0.4, at=seq(5,295, by = 10))
#axis(2, at=seq(0,250, by = 20), las = 2)
box()

temp2 = c(5, 45, 295)

probability = (I0 - koveslighety(temp2,10)) /I0

par(mar = c(4,3,2,2))
for (i in seq(1,3,1)){
 mids = barplot(intensityDistribution(probability[i], I0), 
       main = paste(as.character(temp2[i]), "km"),
       ylab = "probability",
       xlab = "intensity",
       ylim = c(0,1),
       xpd = F,
       xaxs="i", yaxs="i")
  
  axis(1, at=mids, labels = seq(1,9,1)) 
  #axis(lwd=0,side=1,line=-0.4, at=seq(0,1,0.2))
  #axis(2, at=seq(0,250, by = 20), las = 2)
  box()
  #Sys.sleep(10)
}

dev.off()


# Data Generation  Gaussian Distribution in R ####
I0 = 9

numberData = 700

set.seed(12345)
r =  rnorm(numberData, 150, 50)

r[r<0]=0
r[r>300]=300

data = (I0 - koveslighety(r,10))+rnorm(numberData,0,1)

Data = data.frame(Is = round(data), R = r)

Data$Is[Data$Is > I0] = I0
Data$Is[Data$Is < 1] = 1


save(Data, file = "data/syntheticDataNormal.RData")


# Data Generation  lognormal near Distribution in R ####
I0 = 9

numberData = 700

set.seed(12345)
r =  rlnorm(numberData, log(50), log(2))


r[r<0]=0
r[r>300]=300

data = (I0 - koveslighety(r,10))+rnorm(numberData,0,1)



Data = data.frame(Is = round(data), R = r)

Data$Is[Data$Is > I0] = I0
Data$Is[Data$Is < 1] = 1

save(Data, file = "data/syntheticDataLog1.RData")

# Data Generation  lognormal distance Distribution in R ####
I0 = 9

numberData = 700

set.seed(12345)
r =  abs(300 -rlnorm(numberData, log(50), log(2)))


r[r<0]=0
r[r>300]=300

data = (I0 - koveslighety(r,10))+rnorm(numberData,0,1)


Data = data.frame(Is = round(data), R = r)


Data$Is[Data$Is > I0] = I0
Data$Is[Data$Is < 1] = 1

save(Data, file = "data/syntheticDataLog2.RData")

# Figure Koveslighety ####

pdf(file = "../thesis/Figures/koveslighety.pdf",
    height = 9,
    width = 16,
    pointsize = 25)


plot(r, data,
     xlim = c(0,300),
     ylim = c(0,10),     
     main = "synthetic Intensity Data",
     ylab = expression("Intensity at Site I"[" s"]),
     xlab = "Distance R [km]",
     xaxt = "n",
     pch = 19,
     lwd = 5) 

axis(1, at=seq(0,max(Data$R)+1, by = 10), labels = NA)
axis(lwd=0,side=1,line=-0.4, at=seq(0,max(Data$R)+1, by = 20))
#axis(2, at=seq(0,250, by = 20), las = 2)
box()

dev.off()

# Figure synthetic Data ####
pdf(file = "../thesis/Figures/syntheticIntensities.pdf",
    height = 9,
    width = 16,
    pointsize = 25)

plot(Data$R, round(Data$Is),
     xlim = c(0,300),
     ylim = c(1,10),     
     main = "synthetic Intensity Data",
     ylab = expression("Intensity at Site I"[" s"]),
     xlab = "Distance R [km]",
     xaxt = "n",
     pch = 19) 

axis(1, at=seq(0,max(Data$R), by = 10), labels = NA)
axis(lwd=0,side=1,line=-0.4, at=seq(0,max(Data$R)+1, by = 20))
#axis(2, at=seq(0,250, by = 20), las = 2)
box()

dev.off()

# Figure Prior ####

distance = seq(0,300,10)


priorMean = (1/(1+distance/60))^(1/I0)
#priorMean = (I0 -koveslighety(distance,10))/I0
#priorMean = exp(-0.009*distance)
#priormean = 1/300*distance 
maxVariance = priorMean*(1-priorMean)

bound1 = priorMean*(1-priorMean)^2/(2-priorMean)
bound2 = priorMean^2*(1-priorMean)/(1+priorMean)


priorVariance = maxVariance *0.09

#priorVariance = 0.99*bound1
#priorVariance[1]= (bound1[1]+bound2[1])/2

priorVariance = 0.007
alpha =  ((1 - priorMean) / priorVariance - 1 / priorMean) * priorMean ^ 2
beta = alpha * (1 / priorMean - 1)


sdPrior =  sqrt((alpha*beta)/((alpha + beta + 1)*(alpha+beta)^2))

pdf(file = "../thesis/Figures/tsapanosPrior.pdf",
    height = 9,
    width = 16,
    pointsize = 25)

layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow = TRUE))

plot(distance, priorMean, xlim = c(0,300),
     ylim = c(0,1),     
     main = "Prior Distribution",
     ylab = "probability",
     xlab = "Distance R [km]",
     xaxt = "n",
     pch = 19, xaxs="i", yaxs="i")
lines(distance, priorMean, lwd = 5)

segments(distance, priorMean-sdPrior,distance, priorMean+sdPrior)
epsilon = 2
segments(distance-epsilon,priorMean-sdPrior,distance+epsilon,priorMean-sdPrior)
segments(distance-epsilon,priorMean+sdPrior,distance+epsilon,priorMean+sdPrior)

axis(1, at=seq(0,max(Data$R)+1, by = 10), labels = NA)
axis(lwd=0,side=1,line=-0.4, at=seq(0,max(Data$R)+1, by = 20))
#axis(2, at=seq(0,250, by = 20), las = 2)
#box()

temp2 = c(20, 150, 280)



for (i in seq(1,3,1)){
  plot(seq(0,1,0.01), dbeta(seq(0,1,0.01),alpha[temp2[i]/10], beta[temp2[i]/10]), 
       main = paste(as.character(temp2[i]), "km"),
       ylab = "density",
       xlab = "probability p",
       pch = 19,
       xaxt = "n",
       type = "l",
       lwd = 5)
  
  axis(1, at=seq(0,1,0.1), labels = NA) 
  axis(lwd=0,side=1,line=-0.4, at=seq(0,1,0.2))
  #axis(2, at=seq(0,250, by = 20), las = 2)
  #box()
  #Sys.sleep(0.1)
}

dev.off()
