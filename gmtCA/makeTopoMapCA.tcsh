#!/bin/sh

gmt set FORMAT_GEO_MAP ddd:mm:ssF
gmt set FONT_ANNOT_PRIMARY 11p
gmt set FONT_LABEL 14p

region="-R67/80/38/44"
boundaries="-Ba2g2"
projection="-Jm1"
misc="-K -Y10"
post="CentralAsia.ps"

#grdcut topo15.grd -GCentralAsia.grd $region -V
#grdgradient CentralAsia.grd -A0 -Ne0.5 -fg -GCentralAsiaGradient.grd
#grd2cpt CentralAsia.grd -Crelief > CentralAsia.cpt
#makecpt -Ctopo -T-12000/8000/50 -Z > CentralAsia.cpt

gmt psbasemap $region $projection $boundaries $misc > $post

gmt grdimage CentralAsia.grd -ICentralAsiaGradient.grd -R -J -B -CCentralAsia.cpt -K -P -O >> $post

#grdcontour CentralAsia.grd -C500 -R -J -O -V -K >>$post

gmt pscoast -R -J -B -N1/0.75p,0/0/0 -Di -K -O >> $post


cat CACountries.txt | gmt pstext -R -J -O -K >> $post
cat CACapitals.txt | gmt pstext -R -J -Dj0.12c/0.12c -O -K >> $post

cat CAEQ.txt | gmt psxy -J -R -Sa0.5c -G204/0/0 -W0.25p -O -K >> $post
cat CACapitals.txt | gmt psxy -J -R -Ss0.25c -G204/0/0 -W0.25p -O >> $post

#set misc = "-K -O -Bnesw+g255"
#psbasemap -R65/82/35/45 -Jm1 $misc -Y2.0i >> $post

#psscale -CCentralAsia.cpt -D0.5c/2.5c/4c/0.4c -B:"m": -L -O >> CentralAsia.ps

ps2pdf CentralAsia.ps CentralAsia.pdf
pdfcrop CentralAsia.pdf ../thesis/Figures/CentralAsia.pdf
