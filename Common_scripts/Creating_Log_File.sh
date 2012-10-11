#!/bin/bash
# This script is to be run from the Photometry_Same_Field_convolved.py
# 1st argument is image name, 2nd argument is Goodstars coordinate file and 3rd argument is Source.coo and 4th arg is Output filename
# It will create the Log file containg the continious photometry results. (only data after DATAMIN is written by this script.
# Fields upto DATAMIN is written from the Photometry python script.
# Header of the ./Photometry_Log.txt will be
# Image FWHM FILTER EXPTIME TM_START SKY SKYSIGMA DATAMIN Nebulav1647Mag NebulaMag MagsQphotinput.txt.. X Y V1647Mag X Y Mag X Y Mag X Y Mag ...
#-----------------------------indiajoe@gmail.com

imgfile=$1
#echo -n "Magnitude of Nebula including v1647 = " >> ./Photometry_Log.txt
echo -n "`tail -1 $imgfile".mag.2" | gawk '{print $5}'`  " >> $4
#echo -n "Magnitude of Nebula = " >> ./Photometry_Log.txt
for nebulamagfiles in $imgfile"_Intp.fits.mag."*
do 
    echo -n "`tail -1 $nebulamagfiles | gawk '{print $5}'`  " >> $4
done
#echo -n "`tail -1 $imgfile"_Intp.fits.mag.1" | gawk '{print $5}'`  " >> $4
#Now finding the mags of Source.coo contents from the als file.
{ while read line
    do
    Xa=`echo $line | gawk '{print $1}'`
    Ya=`echo $line | gawk '{print $2}'`
    echo -n "$Xa  $Ya  " >>  $4
  # To find the magnitude of Source from the als file 
    echo -n "`gawk '{dx=($2 - xa)**2 ; dy=($3 - ya)**2 ; if ((dx <= 4) && (dy <= 4)) {print $4  }}' xa=$Xa ya=$Ya $imgfile".als.1"`  "  >>  $4
    done } < $3
#Now finding the mags for remaing good stars. Star A is the first entry among them.
{ while read line
    do
        Xa=`echo $line | gawk '{print $1}'`
	Ya=`echo $line | gawk '{print $2}'`
	echo -n "$Xa  $Ya  " >>  $4
    # To find the magnitude of Star from the als file echo 
	echo -n "`gawk '{dx=($2 - xa)**2 ; dy=($3 - ya)**2 ; if ((dx <= 4) && (dy <= 4)) {print $4 }}' xa=$Xa ya=$Ya $imgfile".als.1"`  "  >>  $4

    done } < $2

#Now reading out the phot results of stars in the image optained by subtracting interpolated image. (for McNeil only)
#---- Done from the python photmetry script itself.

# echo >> $4