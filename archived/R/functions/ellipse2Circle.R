
ellipse2Circle = function(F1,F2, x,y){
  
  r1 = sqrt((F1[1]-x)^2 + (F1[2]-y)^2)
  r2 = sqrt((F2[1]-x)^2 + (F2[2]-y)^2)
  c = sqrt(sum((F1-F2)^2))/2
  a = (r1+r2) /2
  b = sqrt(a^2-c^2)

toDegrees = 180/pi

azimuth = asin(abs((F1-F2))[1]/(2*c))*toDegrees 
theta = 90 - azimuth 

x1 = x
y1 = y
x2 = cos(-theta)*x1-sin(-theta)*y1 
y2 = sin(-theta)*x1+cos(-theta)*y1 

x3 = x2 *b/a
y3 = y2

gamma = atan(y3/x3)-atan(y2/x2)+abs(theta)
x4 = cos(gamma)*x3-sin(gamma)*y3
y4 = sin(gamma)*x3+cos(gamma)*y3

return(list(x4,y4))
}