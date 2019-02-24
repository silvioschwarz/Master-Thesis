#------------------------------------------------------#
#
#   Binomial Distribution
#   Compute binomial Distribution of Intensities
#
#   Silvio Schwarz
#   (C) 2014-2015 Silvio Schwarz. GNU GPL 3.
#------------------------------------------------------#



intensityDistribution = function(probability, I0){
  # compute probability distribution
  
  probability[probability > 0.98] = 0.98
  
  distribution = dbinom(x = seq(1, I0, 1), size = I0, prob = probability)
  #modee = which.max(prob)
  #meann = prob * seq(1, I0, 1)
  
  return(distribution)
  #return(append(prob,mode))
}