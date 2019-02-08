#------------------------------------------------------#
#
#   Load Packages
#   load packages and check whether installed
#   if not installed, install.
#
#   Silvio Schwarz
#   (C) 2014-2016 Silvio Schwarz. GNU GPL 3.
#------------------------------------------------------#


loadPackages = function(packageList){
  if (length(setdiff(packageList, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packageList, rownames(installed.packages())))  
  }
  
  sapply(packageList, require, character.only = T)
  
}