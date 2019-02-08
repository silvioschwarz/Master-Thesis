gmt5_loc=/usr/local/lib
gmt5_bin=gmt5
gmt5_dev=gmt5-dev
gmt5_dat=gmt5-data

sudo rm -r $gmt5_loc/$gmt5_bin

sudo apt-get install -y --allow-unauthenticated subversion build-essential cmake libcurl4-gnutls-dev libgdal1-dev libfftw3-dev libpcre3-dev libnetcdf-dev liblapack-dev libblas-dev graphicsmagick texlive texlive-latex-extra python-sphinx

#building

cd ~/Downloads

mkdir $gmt5_bin
cd $gmt5_bin


svn checkout svn://gmtserver.soest.hawaii.edu/gmt/trunk $gmt5_dev

mkdir $gmt5_dat
cd $gmt5_dat

download1=ftp://ftp.soest.hawaii.edu/gshhg/ 
GSHHG=$(wget -q $download1 -O - | egrep -o "gshhg-gmt-[0-9\.]+.tar.gz"| sort -V  | tail -1)
wget -c -O gshhg-gmt.tar.gz $download1$GSHHG

download2=ftp://ftp.soest.hawaii.edu/dcw/
DCW=$(wget -q $download2 -O - | egrep -o "dcw-gmt-[0-9\.]+.tar.gz"| sort -V  | tail -1)
wget -c -O dcw-gmt.tar.gz $download2$DCW

cat *.tar.gz | tar -zxvf - -i

sudo rm -r *.tar.gz
GSHHG=$(ls | tail -1)
DCW=$(ls | head -1)

cd ..

cp $gmt5_dev/cmake/ConfigUserTemplate.cmake $gmt5_dev/cmake/ConfigUser.cmake

## change 3 items in ConfigUser.cmake
echo "   -> modify ConfigUser.cmake" 
 ## 1 - installation directory
sed -i -e "s/.*CMAKE_INSTALL_PREFIX.*prefix_path.*/set (CMAKE_INSTALL_PREFIX \"${gmt5_loc////\\/}\/$gmt5_bin\")/g" $gmt5_dev/cmake/ConfigUser.cmake
 ## 2 - path to gshhg files

sed -i -e "s/.*GSHHG_ROOT.*/set (GSHHG_ROOT \"${gmt5_loc////\\/}\/$gmt5_bin\/$gmt5_dat\/$GSHHG\")/g" $gmt5_dev/cmake/ConfigUser.cmake
 ## 3 - path to dcw files

sed -i -e "s/.*DCW_ROOT.*/set (DCW_ROOT \"${gmt5_loc////\\/}\/$gmt5_bin\/$gmt5_dat\/$DCW\")/g" $gmt5_dev/cmake/ConfigUser.cmake

cd ~/Downloads
sudo mv gmt5/ /usr/local/lib/gmt5

cd $gmt5_loc/$gmt5_bin/$gmt5_dev/
sudo mkdir build
cd build
sudo cmake ..
sudo make

sudo make docs_man
sudo make docs_html
sudo make docs_pdf
sudo make install

## add GMT5/bin to .bashrc PATH
if [ $( grep $gmt5_bin $HOME/.bashrc &> /dev/null; echo $? ) -ne 0 ]; then
  paths+=( $gmt5_loc/$gmt5_bin/bin ); fi
fi

