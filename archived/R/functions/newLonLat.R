newLonLat =
  function (lon, lat, bearing, distance) 
  {
    rad = pi/180
    a1 = lat * rad
    a2 = lon * rad
    tc = bearing * rad
    d = distance/6378.145
    nlat = asin(sin(a1) * cos(d) + cos(a1) * sin(d) * cos(tc))
    dlon = atan2(sin(tc) * sin(d) * cos(a1), cos(d) - sin(a1) * 
                    sin(nlat))
    nlon = ((a2 + dlon + pi)%%(2 * pi)) - pi
    npts = cbind(nlon/rad, nlat/rad)
    return(npts)
  }