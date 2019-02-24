#------------------------------------------------------#
#
#   Intensity Distribution
#   Compute binomial distribution of macroseismic 
#   intensities
#
#   Silvio Schwarz
#   (C) 2014-2016 Silvio Schwarz. GNU GPL 3.
#------------------------------------------------------#



intensityDistribution = function(probability, I0){
  # compute probability distribution
  
  probability[probability > 0.98] = 0.98
  probability[probability < 0] = 0
  
  distribution = dbinom(x = seq(1, I0, 1), size = I0, prob = probability)
  #modee = which.max(distribution)
  #meann = prob * seq(1, I0, 1)
  
 
  return(distribution)
}