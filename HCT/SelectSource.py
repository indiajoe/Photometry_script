#!/usr/bin/env python
#USAGE:  ./SelectSource.py qphot=1,2,3 source=1 stds=1,2,3 filter=75
import sys
Qphot=[]
Source=[]
Stds=[]
Filter=['6','75','76','11']
if len(sys.argv) < 2 : 
    print('USAGE:  ./SelectSource.py qphot=1,2,3 source=1 stds=1,2,3 filter=75')
    exit()
for arg in sys.argv[1:]:
    if arg[0:6]=='qphot=' : Qphot=arg[6:].split(',')
    elif arg[0:7]=='source=' : Source=arg[7:].split(',')
    elif arg[0:5]=='stds=' : Stds=arg[5:].split(',')
    elif arg[0:7]=='filter=' : Filter=arg[7:].split(',')
    else : print('Unknown input :'+arg)
#print('Searching for qphot=',Qphot,'source=',Source,'stds=',Stds, 'Filter',Filter)
StartStd=22
StartSource=18
StartQphot=8
FinalMagIN=open('FinalDiffMag_Output.txt','r')
for iline in FinalMagIN.readlines():
    iline=iline.rstrip()
    iline=iline.split()
    output=[]
    if len(iline) > 4 and iline[0][0] != '#' and iline[2] in Filter :  #Data line,Not a discarded data line,Required filter
        output=[iline[4],iline[0],iline[2]]  # JD Fitsfile Filter 
        gstars=iline[-2].split('_')
        for star in Qphot :
            output.append(iline[StartQphot+(int(star)-1)*2])
        for star in Source :
            output.append(iline[StartSource+(int(star)-1)*4 +2])
        for star in Stds :
            if star in gstars:
                pos=gstars.index(star)
                output.append(iline[StartStd+pos*4 +2])
                output.append(iline[StartStd+(len(gstars)-1)*4+ (StartStd-StartSource)+ pos*4 +2]) #Only for McNeil
            else :
                output.append('INDEF')
                output.append('INDEF')     #Only for McNeil
        #Print only if atleast one Mag is not INDEF
        for mag in output[3:]:
            if mag != 'INDEF' : 
                print(' '.join(output))
                break

FinalMagIN.close()
            
        
