# Initialization####

source("functions/initialization.R")

packageList = c("caret", "beepr","nls2","parallel","foreach","doParallel", "fields", "geosphere")

initialization(packageList)

# Combining Catalogues ####
library(fields)


event = read.csv("data/CAEQ.csv", sep = "\t")

MLH = 0.47* event$K - 1.15
MW =ifelse(MLH <= 6.1, 0.67*MLH + 2.07,0.99*MLH + 0.08)
#E.M. Scordilis


faultLength = 10^(-3.22+ 0.69*MW)


faultPoints1 = destPoint(cbind(event$LON, event$LAT), event$strike1, faultLength*500)
faultPoints2 = destPoint(cbind(event$LON, event$LAT), abs(180 -event$strike1), faultLength*500)



eventData = data.frame(ID = event$idL,
                       year = event$year,
                       month = event$month,
                       day = event$day,
                       LON = event$LON,
                       LAT = event$LAT,
                       strike = event$strike1,
                       Imax = event$Imax,
                       depth = event$depth,
                       MWHARVARD = event$Mw,
                       K = event$K,
                       MLH = MLH,
                       MW = MW,
                       faultLength = faultLength,
                       F11 = faultPoints1[,1],
                       F12 = faultPoints1[,2],
                       F21 = faultPoints2[,1],
                       F22 = faultPoints2[,2]
)
write.table(eventData, file = "data/event.csv", row.names = FALSE, sep = "\t")

dataCentralAsia = read.csv("data/IntensityCentralAsia.csv", sep = "\t")

dataCentralAsia$I0 = numeric(nrow(dataCentralAsia))
dataCentralAsia$F11 = numeric(nrow(dataCentralAsia))
dataCentralAsia$F12 = numeric(nrow(dataCentralAsia))
dataCentralAsia$F21 = numeric(nrow(dataCentralAsia))
dataCentralAsia$F22 = numeric(nrow(dataCentralAsia))
dataCentralAsia$R = numeric(nrow(dataCentralAsia))
dataCentralAsia$MLH = numeric(nrow(dataCentralAsia))
dataCentralAsia$MW = numeric(nrow(dataCentralAsia))
dataCentralAsia$MWHARV = numeric(nrow(dataCentralAsia))
dataCentralAsia$faultLength = numeric(nrow(dataCentralAsia))

dataCentralAsia$R = diag(rdist.earth(matrix(c(dataCentralAsia$site_lon, dataCentralAsia$site_lat), ncol = 2), 
                                     matrix(c(dataCentralAsia$eve_lon, dataCentralAsia$eve_lat), ncol = 2), miles = FALSE, R = 6371))

for(i in 1: nrow(dataCentralAsia)){
  
  index = min(which(dataCentralAsia$event[i] == eventData$ID))
  
  
  
  dataCentralAsia$I0[i] = eventData$Imax[index]
  dataCentralAsia$MWHARV[i] = eventData$MWHARVARD[index]
  
  dataCentralAsia$faultLength[i]  = eventData$faultLength[index]
  dataCentralAsia$MW[i] = eventData$MW[index]
  dataCentralAsia$MLH[i] = eventData$MLH[index] 
  dataCentralAsia$F11[i] = eventData$F11[index]
  dataCentralAsia$F12[i] = eventData$F12[index]
  dataCentralAsia$F21[i] = eventData$F21[index]
  dataCentralAsia$F22[i]= eventData$F22[index]
}

library(sp)
library(rgdal)
long2UTM <- function(long) {
  (floor((long + 180)/6) %% 60) + 1
}
#Function
LongLatToUTM<-function(args){
  x =args[1]
  y=args[2]
  xy <- data.frame(X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  zone = long2UTM(x)
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",42," ellps=WGS84",sep='')))
  return(as.data.frame(res))
}

# Example
x<-c( -94.99729,-94.99726,-94.99457,-94.99458,-94.99729)
y<-c( 29.17112, 29.17107, 29.17273, 29.17278, 29.17112)

sapply(1:length(x), function(a)LongLatToUTM(cbind(x,y)[a,]))
LongLatToUTM(x,y)

siteUTm = sapply(1:length(dataCentralAsia$site_lon), function(a)LongLatToUTM(cbind(dataCentralAsia$site_lon,dataCentralAsia$site_lat)[a,]))

X = unlist(siteUTm[1,])
Y = unlist(siteUTm[2,])

dataCentralAsia$X = X
dataCentralAsia$Y = Y

transCoord = matrix(NA,nrow(dataCentralAsia),2)
for(i in 1:nrow(dataCentralAsia)){
  dataTrans = dataCentralAsia[i,]
  temp = ellipse2Circle(LongLatToUTM(cbind(dataTrans$F11,dataTrans$F12)),
                        LongLatToUTM(cbind(dataTrans$F21,dataTrans$F22)), 
                        dataTrans$X,dataTrans$Y)
  transCoord[i,1] = temp[[1]][[1]]
  transCoord[i,2] = temp[[2]][[1]]
  
}

eventUTM = sapply(1:length(dataCentralAsia$eve_lon), function(a)LongLatToUTM(cbind(dataCentralAsia$eve_lon,dataCentralAsia$eve_lat)[a,]))

eventX = unlist(eventUTM[1,])
eventY = unlist(eventUTM[2,])

RExt = sqrt(rowSums((transCoord  - cbind(eventX,eventY))^2))/1000


transCoord = matrix(NA,nrow(dataCentralAsia),2)
for(i in 1:nrow(dataCentralAsia)){
  dataTrans = dataCentralAsia[i,]
  temp = ellipse2Circle(c(dataTrans$F11,dataTrans$F12),
                        c(dataTrans$F21,dataTrans$F22), 
                        dataTrans$site_lon,dataTrans$site_lat)
  transCoord[i,1] = temp[[1]]
  transCoord[i,2] = temp[[2]]
  
}

RExt = diag(rdist.earth(matrix(c(transCoord[,1], transCoord[,2]), ncol = 2), 
                        matrix(c(dataCentralAsia$eve_lon, dataCentralAsia$eve_lat), ncol = 2), miles = FALSE, R = 6371))


data = data.frame(ID = dataCentralAsia$event, 
                  siteLAT =dataCentralAsia$site_lat,
                  siteLON =dataCentralAsia$site_lon,
                  eventLAT = dataCentralAsia$eve_lat,
                  eventLON = dataCentralAsia$eve_lon,
                  K = dataCentralAsia$K,
                  mag = dataCentralAsia$mag,
                  depth = dataCentralAsia$depth,
                  Is = dataCentralAsia$intensity, 
                  I0 = dataCentralAsia$I0, 
                  R = dataCentralAsia$R,
                  RExt = RExt,
                  F11 =dataCentralAsia$F11,
                  F12 =dataCentralAsia$F12,
                  F21 =dataCentralAsia$F21,
                  F22 =dataCentralAsia$F22,
                  Length= dataCentralAsia$faultLength,
                  MLH = dataCentralAsia$MLH,
                  MW = dataCentralAsia$MW,
                  MWHARV = dataCentralAsia$MWHARV)

write.table(data, file = "data/dataCentralAsia.dat", row.names = FALSE, sep = "\t")



#prepare data####
data = read.csv("data/dataCentralAsia.dat",sep = "\t")
# round data to make them integer 
dataComp = data
dataComp$I0 = round(data$I0)
dataComp$Is = round(data$Is)

save(dataComp, file="data/CaseStudyRound.RData")

# floor data to make them integer 
dataComp = data
dataComp$I0 = floor(data$I0)
dataComp$Is = floor(data$Is)

save(dataComp, file="data/CaseStudyFloor.RData")

# ceil data to make them integer 
dataComp = data
dataComp$I0 = ceiling(data$I0)
dataComp$Is = ceiling(data$Is)

save(dataComp, file="data/CaseStudyCeil.RData")

# leave uncertain data out 

certainIndex = intersect(which(data$I0 == round(data$I0)),which(data$Is == round(data$Is)))
dataComp = data[certainIndex,]

save(dataComp, file="data/CaseStudyCertain.RData")

# inflate catalogue 4x since I0+ Is uncertain = 4 possibilities
# 1 IS + I0 floor
# 2 IS floor I0 ceil
# 3 Is ceil I0 ceil
# 4 is ceil I0 floor


sliceIndex = nrow(data)
dataComp = rbind(data,data,data, data)
dataComp$Is = c(floor(dataComp$Is[1:(2*sliceIndex)]),ceiling(dataComp$Is[(sliceIndex*2+1):(sliceIndex*4)]))
dataComp$I0 = c(floor(dataComp$I0[1:sliceIndex]),ceiling(dataComp$I0[(sliceIndex+1):(sliceIndex*2)]),
                floor(dataComp$I0[(sliceIndex*2+1):(sliceIndex*3)]),ceiling(dataComp$I0[(sliceIndex*3+1):(sliceIndex*4)]))

save(dataComp, file="data/CaseStudyInflated.RData")

# inflate only where there are uncertain data
dataComp = data
ind = which(dataComp$Is!=round(dataComp$Is))
alreadyClean = dataComp[-ind,]
toClean = dataComp[rep(ind, each = 2),]
toClean$Is = toClean$Is + rep(c(0.5,-0.5), length(ind))

dataComp = rbind(alreadyClean,toClean)

ind = which(dataComp$I0!=round(dataComp$I0))
alreadyClean = dataComp[-ind,]
toClean = dataComp[rep(ind, each = 2),]
toClean$I0 = toClean$I0 + rep(c(0.5,-0.5), length(ind))

dataComp = rbind(alreadyClean,toClean)

save(dataComp, file="data/CaseStudyInflatedUncertain.RData")

# randomly floor or ceil
data = read.csv("data/dataCentralAsia.dat", sep = "\t")
dataComp = data
ind = which(dataComp$Is!=round(dataComp$Is))
alreadyClean = dataComp[-ind,]
toClean = dataComp[ind,]
toClean$Is = toClean$Is + sample(c(0.5,-0.5),length(ind), replace = T)

dataComp = rbind(alreadyClean,toClean)

ind = which(dataComp$I0!=round(dataComp$I0))
alreadyClean = dataComp[-ind,]
toClean = dataComp[ind,]
toClean$I0 = toClean$I0 + sample(c(0.5,-0.5),length(ind), replace = T)

dataComp = rbind(alreadyClean,toClean)

save(dataComp, file="data/CaseStudyRandom.RData")
