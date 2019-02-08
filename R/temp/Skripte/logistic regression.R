# Rotondi 2004 Marche-Umbria Sequence

# Set workinf directory

path = getwd()
setwd(path)

rm(list = ls())

### loading data

#needed for distance calculations
library(fields)
#needed for crossvaliadation
library(caret)
#for plotting
library(lattice)


data = read.csv("/home/silvio/Dokumente/Masterarbeit/Masterarbeitdata/Rotondi 2004/DZ-47/DZ-47.txt", sep = "/t", header = TRUE)

#data$Is = floor(data$Is)
#data$Ix = floor(data$Ix)

# computing delta I = Ix - Is
data$delta_i = data$Ix - data$Is

# computing intensities normalized by  Ix
data$Isnorm = data$Is / data$Ix

# computing the distance from the epicentre 
distance = diag(rdist.earth(matrix(c(data$LON, data$LAT), ncol = 2), 
                            matrix(c(data$LON_epi, data$LAT_epi), ncol = 2), miles = FALSE, R = 6371))
data$distance = distance

color = floor(data$Is/10)+1

plot(data$distance, data$Ix, col = color)

data = data[data$Ix == 11,]

plot(data$distance, data$Is)
data$label = 0
a = 1:11
for (i in 3:10){
data[data$Is <i,]$label = 0
data[data$Is >=i,]$label = 1

plot(data$distance, data$label, xlim = c(0, 500), ylim = c(0, 1))

fit= glm(data$label ~ data$distance, binomial)


coeff = coef(fit)

a[i] = - coeff[[1]]/coeff[[2]]
b = exp(coeff[[1]]+coeff[[2]]*a)/(1+exp(coeff[[1]]+coeff[[2]]*a))

x=seq(0,500,0.1)
plot(x,(exp(coeff[[1]]+coeff[[2]]*x)/(1+exp(coeff[[1]]+coeff[[2]]*x))), col = 'red',xlim = c(0, 500), ylim = c(0, 1))
#plot(x,(exp(coeff[[1]]+coeff[[2]]*x)/(exp(coeff[[1]]+coeff[[2]]*x)+1)^2), col = 'red',xlim = c(0, 500), ylim = c(0, 1))
points(data$distance, data$label, col = 'blue')
title(paste('threshold: ',as.character(i)))
readline(prompt = "")
}





library(MethComp)
# 'True' values 
M = runif(100,0,5)
# Measurements:
x =         M + rnorm(100)
y = 2 + 3 * M + rnorm(100,sd=2)
# Deming regression with equal variances, variance ratio 2.
Deming(x,y)
Deming(x,y,vr=2)
Deming(x,y,boot=TRUE)
bb = Deming(x,y,boot=TRUE,keep.boot=TRUE)
str(bb)
# Plot data with the two classical regression lines
plot(x,y)
abline(lm(y~x))
ir = coef(lm(x~y))
abline(-ir[1]/ir[2],1/ir[2])
abline(Deming(x,y,sdr=2)[1:2],col="red")
abline(Deming(x,y,sdr=10)[1:2],col="blue")
# Comparing classical regression and "Deming extreme"
summary(lm(y~x))
Deming(x,y,vr=1000000)


