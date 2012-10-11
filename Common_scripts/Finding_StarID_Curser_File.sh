#!/bin/bash
#Run by Photometry script after phot with image file as 1st argument and Goodstars coordinate file as 2nd argument.
#This script will find the Star ID from the .mag.1 file and then create the icommands required to run psf non-interactively.
#---------------------------------indiajoe@gmail.com
rm icommands.in 2> /dev/null

imgfile=$1
#IMG=${imgfile%%.fits}
#sed -n /^$IMG/p $imgfile.mag.1 > imgID_temp.sed
sed -n "s/^${imgfile:0:23}/$imgfile /p"  $imgfile.mag.1 > imgID_temp.sed
{ while read line
    do
        Xa=`echo $line | gawk '{print $1}'`
        Ya=`echo $line | gawk '{print $2}'`
        StarID=`gawk '{dx=($2 - xa)**2 ; dy=($3 - ya)**2 ; if ((dx <= 4) && (dy <= 4)) {print $4 }}' xa=$Xa ya=$Ya imgID_temp.sed `
	echo ":a $StarID" >> icommands.in
    done } < $2

echo "f" >> icommands.in
echo "w" >> icommands.in
echo "q" >> icommands.in