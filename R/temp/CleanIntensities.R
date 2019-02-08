#------------------------------------------------------#
#
#   Clean Intensities
#   Turn List of IDP's to integers
#
#   Silvio Schwarz
#   (C) 2014-2015 Silvio Schwarz. GNU GPL 3.
#------------------------------------------------------#

cleanIntensities = function(){
  
  filename_in = file.choose()
  data = read.csv(filename_in,sep = "/t")
  
  ind = which(data$Is!=round(data$Is))
  alreadyClean = data[-ind,]
  toClean = data[rep(ind, each = 2),]
  toClean$Is = toClean$Is + rep(c(0.5,-0.5), length(ind))
  
  data = rbind(alreadyClean,toClean)
  
  
  ind = which(data$I0!=round(data$I0))
  
  if(length(ind)!= 0) {
    alreadyClean = data[-ind,]
    toClean = data[rep(ind, each = 2),]
    toClean$I0 = toClean$I0 + rep(c(0.5,-0.5), length(ind))
    data = rbind(alreadyClean,toClean)
  }
  
  rownames(data) = NULL
  
  file_out = file.choose()
  
  write.table(data, file = file_out, row.names = FALSE)
  
  return(data)
  
}

#Unit Test####
#filename_in = "./data/DZ47.dat" # UnitTest
