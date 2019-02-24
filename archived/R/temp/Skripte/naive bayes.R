
# Loading data and initilization ####

# Set working directory

path = getwd()
setwd(path)

rm(list = ls())

#needed for distance calculations
library(fields)
#needed for crossvalidation
library(caret)
#for plotting
library(lattice)


data = read.csv("./Masterarbeit/Masterarbeitdata/Rotondi 2004/DZ-47.dat", sep = " ", header = TRUE)
#data = read.csv("./Masterarbeit/Masterarbeitdata/Rotondi 2004/DZ-47/DZ-47.txt", sep = "/t", header = TRUE)
#computing the distance from the epicentre 
#R = diag(rdist.earth(matrix(c(data$LON, data$LAT), ncol = 2), 
#                      matrix(c(data$LON_epi, data$LAT_epi), ncol = 2), miles = FALSE, R = 6371))

#data$Is = floor(data$Is)
#data$I0 = floor(data$I0)



# computing delta I = I0 - Is
data$delta_i = data$I0 - data$Is

# computing intensities normalized by  I0
data$Isnorm = data$Is / data$I0

# binning distance
binwidth = 10
data$R = floor(data$R / binwidth) * binwidth + binwidth#/2 

# calculating the maximum and minimum distance in the dataset
maxdist = max(data$R)
mindist = min(data$R)
# number of bins given maximum distance and bin width (!should be integer!)
numbin = floor((maxdist - mindist) / binwidth)

# set distance equal to middle of bins
d = seq(mindist, numbin * binwidth, binwidth)
# define edges of the bins
edges = seq(mindist, maxdist, binwidth)

# set number of fold to be used in crossvalidation
folds = 10

# initialize variables
pj = array(0,dim=c(7,numbin))
pj_new = array(0,dim=c(7,numbin))
var_new = array(0,dim=c(7,numbin))

coeff1 = array(0,dim=c(7,2))
coeff2 = array(0,dim=c(7,2))

abs_error = matrix(0, nrow = 7, ncol = 1)
fiterror1 =  matrix(0, nrow = 7, ncol = 1)
fiterror2 = matrix(0, nrow = 7, ncol = 1)


#exclude eq 45 27/10/1914 because of anomalous high values at long distances
data = data[!data$ID == 45,]


# distribution of I_s over distance # would be useful for bayesian classifiers! ####

for (e in 5:11){
        
        ioall = which(data$I0 %in% c(e-0.5,e,e+0.5))
        dataplot = data[ioall,]
        
        #Sys.sleep(3)
        
        for(i in 2:e){
                
                data_is = dataplot[which(dataplot$Is %in% c(i-0.5,i,i+0.5)),]
                
                
                par(mfrow=c(1,1)) 
                hist(data_is$R,  breaks = edges, freq = FALSE, xlim = c(0, 500), 
                     main = paste('I0:', as.character(e), 
                                  paste(',    Is:', as.character(i),
                                        paste(',    mean:', as.character(mean(data_is$R)),
                                              paste(',    sd:', as.character(sd(data_is$R)))))))
                lines(seq(0,500,0.1), dlnorm(seq(0,500,0.1), mean(log(data_is$R)),sd(log(data_is$R))),
                      col = 'red')
                #plot(density(data_is$R), xlim = c(0, 500), main = '')
                readline(prompt = "")
                #Sys.sleep(0.5)
        }
        
} 

# prior, mean = distance of köveslighty or other preliminary attenuation relationship, sd = ? decreasing from high to low IS, isoseismals???