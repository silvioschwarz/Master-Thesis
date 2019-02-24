# Rotondi 2004 Marche-Umbria Sequence

# Loading data and initilization ####

# Set working directory

path = getwd()
setwd(path)

# Clear the memory 
rm(list = ls())

# read in additional functions
pathnames = list.files(pattern="[.]R$", path="functions", full.names=TRUE)
sapply(pathnames, FUN=source)

packages = c("fields", "xtable")

loadPackages(packages)


# setting "stylesheet"
par(pch = 19, ann = T, cex = 1)

# reading in data
dataPrior = read.csv("./data/DZ47.dat", sep = "t")
dataPost = read.csv("./data/D47.dat", sep = "t")
table(cut(data$I0, seq(4.1, 11.6, 0.5)))



# binning distance
binwidth = 10

# calculating the maximum and minimum distance in the dataset
mindist = floor(min(data$R))
maxdist = 250
data = data[data$R <= maxdist,]

# number of bins given maximum distance and bin width (!should be integer!)
numbin = maxdist / binwidth

# set distance equal to middle of bins
d = seq(mindist + binwidth / 2, numbin * binwidth, binwidth)
# define edges of the bins
edges = seq(mindist, maxdist, binwidth)

#set R to bin max
data$R = floor(data$R / binwidth) * binwidth + binwidth

# Visualization ####

hist(data$Ix, breaks = seq(1:12))
hist(data$I0, breaks = seq(1:12))
hist(data$Is, breaks = seq(1:12))
hist(data$R, breaks = seq(0,ceiling(maxdist),1))

plot(data$R, data$I0, 
     col = adjustcolor(data$Is, alpha.f = 0.3), 
     xlim = c(0, 580))

table(cut(data$I0, seq(4.1, 11.6, 0.5)))
table(cut(data$Is, seq(0.1, 12.1, 0.5)))
table(cut(data$R, seq(0, maxdist, binwidth)))

# development of binomial distribution over distance ----

datamean =  matrix(0, nrow = maxdist/binwidth, ncol = 7)
datasd =  matrix(0, nrow = maxdist/binwidth, ncol = 7)
n = matrix(0, nrow = maxdist/binwidth, ncol = 7)

for (e in 5:11){
  
  x = seq(10, numbin * binwidth, binwidth)
  
  y = seq(0.5, 11.5, 1)
  
  ioall = which(data$I0 %in% c(e-0.5,e,e+0.5))
  dataplot = data[ioall,]
  
  #Sys.sleep(3)
  
  for(dist in x){
    
    datadist = dataplot[dataplot$R == dist,]            
    datamean[(dist-binwidth)/binwidth +1, e-4] = mean(datadist$Is/datadist$I0)
    datasd[(dist-binwidth)/binwidth +1, e-4] = sd(datadist$Is/datadist$I0)
    
    n[(dist-binwidth)/binwidth +1, e-4] = length(datadist$ID)
  }
  
  par(mfrow=c(3,1))
  barplot(n[,e-4], axes = T, ylim = c(0, max(n[,e-4])), space = 0)
  title(paste('epicentral Intensity I0: ',as.character(e)))
  par(mar=c(0,3,1,1))
  plot(x,datamean[,e-4])
  plot(x,datasd[,e-4])
  readline(prompt = "")
  #Sys.sleep(0.5)
}

# relative frequencies for each I0 RAINBOW Image ####
# TODO add transparency

# ID = 58 seems anomalous
#data = data[!data$ID == 58,]

nobs = matrix(0, nrow = 7, ncol = numbin)
for (e in 5:11){
  
  ioall = which(data$I0 %in% c(e-0.5,e,e+0.5))
  dataplot = data[ioall,]
  
  x = seq(0, numbin * binwidth, binwidth)
  
  y = seq(0.5, 11.5, 1)
  
  #freqs = prop.table(table(cut(dataplot$R, x), cut(dataplot$Isnorm, y)),1)
  nums = table(cut(dataplot$R, x), cut(dataplot$Is, y))
  binnums = rowSums(nums)
  nobs[e-4,] = binnums
  #nums/rowSums(nums)[row(nums)] 
  transp = matrix(rep(binnums/sum(binnums),11), nrow = numbin)
  freqs = prop.table(nums,1)
  freqs[is.na(freqs)] = 0
  par(mfrow=c(1,1)) 
  
  image.plot(freqs, axes = FALSE)
  image(transp, axes = FALSE, col = gray.colors(500, start = 1, end = 0, alpha=0.3), add=T)
  
  
  title(paste('epicentral Intensity I0: ',as.character(e)))
  axis(1, at=seq(0,1,length.out = 11), labels = seq(0,500,50))
  axis(2, at=seq(0,1,length.out = 11), labels = seq(1,11,1), las = 2)
  box()
  
  #readline(prompt = "")
  Sys.sleep(5)
}

# mean and sd over distance for each I0 #####
for (e in 1:7){
  par(mfrow=c(2,1))     
  
  d = seq(mindist - binwidth , numbin * binwidth, binwidth)
  
  to45km = seq(-0.47, -6.47, -1) / -0.055
  from45km = seq(1.015, -4.985, -1)/-0.022
  maxdists = c(to45km[to45km<= 45], from45km[from45km > 45])
  
  datameans = datamean[d <= maxdists[e],e]
  datasds = datasd[d <= maxdists[e],e]
  
  de = d[d <= maxdists[e]]
  plot(d, datamean[,e]/(e+4),
       xlab = 'R [km]', 
       ylab = 'mean of Is',
       main = paste('mean'), 
       ylim = c(0, 1),
       xlim = c (0, 240),
       col = 'red', 
       pch = 19)
  #points(de,datameans, pch = 19)
  title(paste('epicentral Intensity I0: ',as.character(e+4)), outer = T)
  
  #         fit = nls(datameans ~ exp(b*de + c) , 
  #                data = data.frame(de, datameans), 
  #                start = c(b = 0.01, c = 4))
  #         fitall = nls(datamean[,e] ~ exp(b*d + c) , 
  #                   data = data.frame(d, datamean[,e]), 
  #                   start = c(b = 0.01, c = 4))
  #         
  x = seq(0, numbin * binwidth,1)
  #lines(x, predict(fit, list(de = x)))
  #lines(x, predict(fitall, list(d = x)), col = 'red')
  #abline(v=maxdists[e],col=1,lty=1)
  
  
  
  plot( d, datasd[,e]/(e+4), 
        xlab = 'R [km]', 
        ylab = 'standard deviation of Is', 
        main = paste('standarddeviation'), 
        ylim = c(0,0.5 ),
        xlim = c (0, 240),
        col = 'red',
        pch = 19)
  #points(de,datasds,pch = 19)
  
  #abline(lm( datasds ~ de, data = data.frame(de, datasds) ))
  #abline(lm( datasd[,e] ~ d, data = data.frame(d, datasd[,e]) ), col = 'red')
  #abline(v=maxdists[e],col=1,lty=1)   
  
  readline(prompt = "")
}

# proof that Rotondi et al don't use middle of the bins #####
ptest = c(6.2, 5.6, 4.7, 4.1)/7
dist = c(15,25,35,45)

prior5 = c(212, 191, 166, 170,183, 174, 168,175,117,129,143,112,131,157)/223
d5 = c(10,20,30,40,50,60,70,80,90,100,110,130,190,210)

prior6 = c(202,162,134,124,94)/223
d6 = seq(10,50,10)

prior7 = c(205,172,162,141,113)/223
d7 = seq(10,50,10)

prior8 = c(200,177,117,143)/223
d8 = seq(10,40,10)

prior9 = c(196,177,148,129)/223
d9 = seq(10,40,10)

prior1011 = c(195,179,161,143,143)/223
d1011 = seq(10,50,10)

coeff_paper = matrix(c(7.116, 8.207, 8.317, 7.372, 7.026, 5.903, 5.903, 
                       0.160, 0.412, 0.305,  0.318, 0.284, 0.206, 0.206),nrow = 7, ncol = 2)

x = 0:250
plot(d1011, prior1011, xlim = c(0, 250), ylim = c(0,1))
fittest = nls(prior1011 ~ (g1/d1011)^g2, data = data.frame(cbind(d1011, prior1011)), start = list(g1 = 4, g2 = 0.3))
lines(x, predict(fittest, data.frame(d1011=x)))
fittest

myfit = matrix(c(7.4869, 8.2268, 8.1355, 7.3532, 6.8743, 5.6419, 
                 0.1611, 0.4077, 0.3057,  0.3183, 0.2796, 0.2045),nrow = 7, ncol = 2)
coeff_paper = matrix(c(7.116, 8.207, 8.317, 7.372, 7.026, 5.903, 5.903, 
                       0.160, 0.412, 0.305,  0.318, 0.284, 0.206, 0.206),nrow = 7, ncol = 2)


# alternative prior computation ----

test = table(cut(data$I0, seq(4.6, 11.6, 0.5),right=F), 
             cut(data$Is, seq(0.6,12,0.5),right=F)), cut(data$R, seq(0,240,binwidth)))

testarray = array(0, dim = c(7, numbin))

for(i0 in seq(1,7,1)){
  
  temp = 0.5 * test[2*i0-1,,] + test[2*i0 ,,] + 0.5 * test[2*i0+1 ,,]
  temp2 = 0.5 * temp[,2 * (i0 + 4) - 1] + temp[,2 * (i0 + 4)] + 0.5 * temp[,2 * (i0 + 4) + 1]
  temp3 = 0.5 * sum( test[2*i0-1,,]) + sum( test[2*i0 ,,]) + 0.5 * sum( test[2*i0 + 1,,])
  testarray[i0,] =  (temp2/temp3) ^ (1 / (i0 + 4))
  
}



# evolution of beta distribution ####

meanalphaj = colMeans(alphaj[5,,])
meanbetaj = colMeans(betaj[5,,])
for (i in 1:numbin){
  x = seq(0,1, 0.01)
  plot(x, dbeta(x, alphaj0[5,i], betaj0[5,i]), type = "l", ylim = c(0, 30))
  
  lines(x, dbeta(x, meanalphaj[i], meanbetaj[i]), col = "red")
  
  title(paste('epicentral distance: ',as.character(i*10)))
  readline(prompt = "")
}


meanfiterror3 = rowMeans(fiterror3)


# Distributions ####
# Binomial Distribution

binom_distro = function(R,I0){
  # compute probability distribution
  probability = (coeff_post[I0-4,1]/R)^coeff_post[I0-4,2]
  
  probability[probability > 0.98] = 0.98
  
  prob = dbinom(x = seq(1, I0, 1), size = I0, prob = probability)
  #prob_is = prob[floor(Is)]
  mode = which.max(prob)
  
  return(append(prob,mode))
  
}


for(i in 6:11){
  validation_data = which(data.post$I0 %in% c(i-0.5,i,i+0.5))
  data_val = data.post[validation_data,]
  I0 = i
  test = sapply(data_val$R, binom_distro, I0 = I0)
  discr = floor(data_val$Is) - test[dim(test)[1],]
}



prob_bin = array(0,dim=c(7,numbin,12))

for (i in 5:11){
  for (d in seq(5,500,binwidth)){
    
    prob_bin[i-4,(d+5)/binwidth,1:i] = binom_distro(d,i)[1:i]
    
    barplot(prob_bin[i-4,(d+5)/binwidth,1:i], 
            names.arg=seq(1,i,1), 
            xlim = c(0,11.5), ylim = c(0,1))
    title(paste('epicentral distance:', as.character(d), paste(',    I0:', as.character(i))))
    readline(prompt = "")
  }
  
}


io = 6
image.plot(prob_bin[io-4,,], axes = FALSE)     
title(paste('epicentral Intensity: ',as.character(io)))
axis(1, at = seq(0,1,length.out = 11), labels = seq(0,500,50));
axis(2, at = seq(0,1,length.out = 12), labels = seq(1,12,1), las = 2);
box()

# Predictive Distribution

pred_distro = function(IO, IS, alpha, beta){
  prob = choose(I0,IS) * 
    (gamma(alpha + beta) / (gamma(alpha) * gamma(beta))) *
    ((gamma(alpha + IS)*gamma(beta + I0 - IS)) / gamma(alpha + beta + I0))
  
}

prob_pred = array(0, dim = c(7,numbin,12))

for (i0 in 5:11){
  for (i in 1:i0){
    for (j in seq(1,numbin,1)){
      
      alpha = mean(alphaj[i0-4,,j])
      beta = mean(betaj[i0-4,,j])
      
      prob_pred[i0 - 4,j,i] = pred_distro(i0,i,alpha,beta)
    }
  }
}


io = 9
image.plot(prob_pred[io-4,,], axes = FALSE)     
title(paste('epicentral Intensity: ',as.character(io)))
axis(1, at=seq(0,1,length.out = 11), labels = seq(0,500,50));
axis(2, at=seq(0,1,length.out = 12), labels = seq(1,12,1), las = 2);
box()


# Plotting the distributions 

for (i in 5:11){
  par(mfrow=c(1,1))  
  
  Sys.sleep(5)
  for (d in seq(5,500,binwidth)){
    
    #result[i-4,(d+5)/binwidth,1:i] = dbinom(x = seq(1, i, 1), size = i, prob = (coeff_post[i-4,1]/d)^coeff_post[i-4,2])
    
    
    
    barplot(prob_bin[i-4,(d+5)/binwidth,], 
            names.arg=c(1,2,3,4,5,6,7,8,9,10,11,12), 
            xlim = c(0,11.5), ylim = c(0,1))
    title(paste('epicentral distance:', as.character(d), paste(',    I0:', as.character(i))))
    #readline(prompt = "")
    Sys.sleep(1)
    
  }
  #         image.plot(result[i-4,,], axes = FALSE)     
  #         
  #         
  #         title(paste('epicentral Intensity: ',as.character(i)))
  #         axis(1, at=seq(0,1,length.out = 11), labels = seq(0,500,50));
  #         axis(2, at=seq(0,1,length.out = 11), labels = seq(1,11,1), las = 2);
  #         box()
  #         
  #         readline(prompt = "")
  
}

# Validation ####

binom_distro = function(R,I0){
  # compute probability distribution
  probability = (coeff_paper_post[I0-4,1]/R)^coeff_paper_post[I0-4,2]
  
  probability[probability > 0.98] = 0.98
  
  prob = dbinom(x = seq(1, I0, 1), size = I0, prob = probability)
  #prob_is = prob[floor(Is)]
  mode = which.max(prob)
  
  return(append(prob,mode))
  
}


pred_distro = function(R, I0, IS){
  alpha = alphaj[I0-4, floor(R/binwidth)+1]
  beta = betaj[I0-4, floor(R/binwidth)+1]
  
  
  prob = choose(floor(I0),floor(IS)) * 
    (gamma(alpha + beta) / (gamma(alpha) * gamma(beta))) *
    ((gamma(alpha + IS)*gamma(beta + I0 - IS)) / gamma(alpha + beta + I0))
  
  return(prob)
}


# validation data  = camerino earthquake 28/7/1799 I0 = 9

data.post = read.csv("./Masterarbeit/Masterarbeitdata/Rotondi 2004/D47.dat", sep = "t", header = TRUE)
data_val = data.post[data.post$ID ==  9,]


plot(data_val$R, data_val$Is)


# plug-in binomial distribution
bin_prob = sapply(data_val$R,binom_distro, I0 = data_val$I0[1])
mode_bin = bin_prob[nrow(bin_prob),]
prob_bin_is = diag(bin_prob[floor(data_val$Is),])
prob_bin_ipred = diag(bin_prob[mode_bin,])

# predictive distribution
pred_prob = pred_distro(data_val$R, data_val$I0, data_val$Is)


mode_pred = pred_prob[nrow(bin_prob),]
prob_pred_is = diag(pred_prob[floor(data_val$Is),])
prob_pred_ipred = diag(pred_prob[mode_pred,])

# logarithmic scoring rule
# -1/N log prod(posterior Probs)

scoring_bin = -log(prod(prob_bin_is))/length(prob_bin_is)


# Odds
# Prob(Is)/Prob(Ipred)
odds_bin = -log(prod(prob_bin_is/prob_bin_ipred))/length(prob_bin_ipred)
# absolute descrepancy
# 1/N abs(Is - Ipred)

discr_bin = sum(abs(data_val$Is - mode_bin))/length(mode_bin)

# precision, accuracy, sensitivity, specificity...
#      Ipred
# Is   1 2 3 4 5 6 7 8 9 10 11 12
#      1 9 2 6 3 . . . . . .  .   .   specificity = # Ipred = i / #Is = i ??? diogonal / sum of all??
#      2 8 4 2 7
#      3 7 2 8 5
#      4 6 6 3 7
#      5 5 2 6 2
#      6 4 6 2 6
#      7 3 1 7 2
#
#
#
#
#
#



### beta binomial distro = optimal bayes classifier = posterio predictive distro = averaging over posterior p * binomial likelihood?
# = integral likelihood * posterior d parameter
# binomial likeliehood p(D|p) = (n over k) p^k * (1-p)^n-k
# beta prior p(p) = p^alpah-1*(1-p)^beta-1 normalized by betafunction(alpha,beta)
# posterior = p(p|D) = p(D|p)*p(p) = beta posterior = (n over K) / betafunctio(alpha,beta) * p^k+alpha-1 * (1-p1)^n-k+beta-1
# n over k = gammaf(n +1) /(gammaf(k +1) gammaf(n-k+1))
# betafunction(alpha,beta) = gammaf(alpha) + gammaf(beta)/gammaf(alpha +beta)



