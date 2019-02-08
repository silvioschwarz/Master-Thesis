# Clear the memory 
rm(list = ls())

# Set working directory

path = getwd()
setwd(path)


# needed for distance calculations
library(fields)

d47 = "./Masterarbeit/Masterarbeitdata/Rotondi 2004/D47.csv"
dz47 = "./Masterarbeit/Masterarbeitdata/Rotondi 2004/DZ47.csv"

data = read.csv(dz47)
data_post = read.csv(d47)

# computing the distance from the epicentre 
data$R = diag(rdist.earth(matrix(c(data$LON, data$LAT), ncol = 2), 
                     matrix(c(data$LON_epi, data$LAT_epi), ncol = 2), miles = FALSE, R = 6371))

data_post$R = diag(rdist.earth(matrix(c(data_post$LON, data_post$LAT), ncol = 2), 
                          matrix(c(data_post$LON_epi, data_post$LAT_epi), ncol = 2), miles = FALSE, R = 6371))

# merging the date
date = character(length(data$R))

for(i in seq(1, length(data$R), 1)){
date[i] = paste(c(data$Da[i],data$Mo[i],data$Ye[i]), collapse = ".")
}
data$Date = date

test = data[c("ID","Date","LAT","LON","LAT_epi","LON_epi","R","Ix","I0","Is")]

write.table(test, file = "./Masterarbeit/Masterarbeitdata/Rotondi 2004/DZ47.dat", row.names = FALSE, sep = "/t")

# posterior data
date.post = character(length(data_post$R))

for(i in seq(1, length(data_post$R), 1)){
        date.post[i] = paste(c(data_post$Da[i],data_post$Mo[i],data_post$Ye[i]), collapse = ".")
}
data_post$Date = date.post

test2 = data_post[c("ID","Date","LAT","LON","LAT_epi","LON_epi","R","Ix","I0","Is")]
write.table(test2, file = "./Masterarbeit/Masterarbeitdata/Rotondi 2004/D47.dat", row.names = FALSE, sep = "/t")








# compute delta I
data$delta_i = data$I0 - data$Is


plus = data[data$I0 != round(data$I0),]
plus$I0.new = plus$I0 + 0.5
minus = data[data$I0 != round(data$I0),]
minus$I0.new = minus$I0 - 0.5
zero = data[data$I0 == round(data$I0),]
zero$I0.new = zero$I0 

data = rbind(zero, plus, minus)

plus = data[data$Is != round(data$Is),]
plus$Is.new = plus$Is + 0.5
minus = data[data$Is != round(data$Is),]
minus$Is.new = minus$Is - 0.5
zero = data[data$Is == round(data$Is),]
zero$Is.new = zero$Is 

data = rbind(zero, plus, minus)


plot(data$R, data$I0, col = data$Is)







