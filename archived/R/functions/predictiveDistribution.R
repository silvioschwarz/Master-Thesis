
predictiveDistribution = function(alphaj,betaj,i0){
  i = seq(0,i0,1)
  returnval =choose(i0,i)*gamma(alphaj+betaj)/(gamma(alphaj)+gamma(betaj))* 
    gamma(alphaj+i)*gamma(betaj+i0-i)/gamma(alphaj+betaj+i0)
  
  returnval = returnval/sum(returnval)
  returnvalue = returnval[2:(i0+1)] +rep(returnval[1]/i0,i0)
  return(returnvalue)
}
