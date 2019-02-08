#------------------------------------------------------#
#
#   Köveslighety Intensity Attenuation
#   Computes intensity attenuation following Köverslighety
#
#   Silvio Schwarz
#   (C) 2014-2015 Silvio Schwarz. GNU GPL 3.
#------------------------------------------------------#

koveslighety = function(r,h){
  
  dh = sqrt(r^2+h^2)
  
  data = 3*log10(dh/h) + 3 * 0.002 *log(dh-h)
  
  return(data)
}