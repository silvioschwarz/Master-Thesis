# Initialization####

source("functions/initialization.R")

packageList = c("caret", "beepr","nls2","parallel","foreach","doParallel", "xtable")

initialization(packageList)

# event Data####
event = read.csv("data/event.dat", sep = "\t")
xtable(event[,-c(1,10)],digits=c(0,0,0,0,2,2,1,0,0,1,3,2,2,3,3,3,3))

# visualizations ####

data = read.csv("data/dataCentralAsia.dat", sep = "\t")

hist(data$Is, breaks = seq(0.5,9.5,1))
hist(data$R)
hist(data$I0, breaks = seq(4.5,9.5,1))

# How much uncertain Data?
# I_s
length(which(data$Is != round(data$Is))) / length(data$Is) * 100

# I_max
length(which(data$I0 != round(data$I0))) / length(data$I0) * 100

# both uncertain
length(which(data$I0 != round(data$I0) & data$Is != round(data$Is))) / length(data$I0) * 100

# Tables ####

loadPackages("xtable")
tab = apply(table(cut(data$I0, seq(9.5,5,-0.5), right = T),
                  cut(data$Is, seq(1.5,9.5,0.5),right = T)),2,rev)
tab = rbind(tab, colSums(tab))
tab = cbind(tab, rowSums(tab))
rownames(tab) = c(seq(9,5,-0.5),"total")
colnames(tab) = c(seq(1.5,9,0.5),"total")
xtable(tab, digits=rep(0, ncol(tab)+1))
