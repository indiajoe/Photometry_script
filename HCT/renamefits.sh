#!/bin/sh
ls */* -1 > images
for img in `cat images`
do
   mv $img $img`echo ".fits"`
done

