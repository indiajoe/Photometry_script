#!/bin/sh
# This script is to suffix the .fits extension to the extension less hct raw files
ls */* -1 > images
for img in `cat images`
do
   mv $img $img`echo ".fits"`
done

