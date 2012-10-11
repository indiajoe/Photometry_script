#!/usr/bin/python
# IMP:  Keep ds9 open 
# This script will guide you through the V R and I images directories.
# Tell the script which combination of images to align , combine etc.
# After that you can do photometry.
# Enjoy!!!---------------------- indiajoe@gmail.com

import os
import glob
from pyraf import iraf
import pyfits
import sys, traceback 
dirs=open('directories','r')
for night in dirs.readlines() :
    night=night.rstrip()
    iraf.cd(night) 
    print night
    filters=glob.glob('image[RVI]')    #Only going to image[RVI] 
    filters.sort()
    
    for imgdir in filters :
        iraf.cd(imgdir)
        while 1 :
            print night+'   '+imgdir
            images=glob.glob('[sn][nz]*.fits')
            images.sort()
            imgs_txt=open('imgs.txt','w')
            i=0
            for img in images:
                hdulist=pyfits.open(img)
                Exptime=hdulist[0].header.get('EXPTIME')
                Commentx=hdulist[0].header.get('COMMENTX')
                Obj=hdulist[0].header.get('OBJECT')
                Comment=hdulist[0].header.get('COMMENT')
                Filter=hdulist[0].header.get('FILTER')
                hdulist.close()
                imgs_txt.write(str(i)+' '+img+' '+str(Obj)+' '+str(Exptime)+' '+str(Commentx)+' '+str(Comment)+' '+str(Filter)+'\n')
                i=i+1
        #   Now list is ready, continuing with what to do
            imgs_txt.close()
#            os.system('../../Display_img_list.sh  &')    # Uses kdialog to display the text file
            os.system('zenity --text-info --title "Images in '+night+'/'+imgdir+'" --filename=imgs.txt --width 600 --height 250 &')

            print ('What do you want to do? Enter a = Align ; c = Combine ; i = Inspect ; e = Exit and go to Next filter directory')
            choice=raw_input('Enter your Choice (a,c,i or e) :')
            if ( choice == "e" ) : break   # breaking to move to next filter image directory
            elif (choice == "i") :
                while 1 :
                    imgEX=input('Enter the Sl.No# of Image you want to Examine [Enter a large number not in image list to exit] : ')
                    if imgEX >= len(images) : break
                    iraf.display(images[imgEX],1)
                    iraf.imexam()
            elif (choice == "c" ) :  # Simply combining 
                inpVAR=""
                while 1 :
                    imgCOM=input('Enter the Sl.No# of next Image to combine [Enter a large number not in image list to exit] : ')
                    if imgCOM >= len(images) : break
                    inpVAR=inpVAR+','+images[imgCOM]
                inpVAR=inpVAR[1:]
                name=raw_input('Enter the name for combined image without .fits extension : ')
                iraf.imcombine(input=inpVAR, output=name+'.fits',combine="average")
            elif ( choice == "a" ) :    # Moving forward to do alignment
                inp=input('Enter the Sl.No# of Reference Image :')
                Refimage=images[inp]
                iraf.display(Refimage,1) 
                print ('Press _a_ over some good stars to align, u can press s also, but DONT press r \n')
                imx=iraf.imexam(Stdout=1)
                Xref=eval(imx[2].split()[0])    #Using the first star as ref star for cose shift calculation
                Yref=eval(imx[2].split()[1])
                foo=open(Refimage+'.coo','w')
                i=2
                while i < len(imx) :               
                    foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
                    i=i+2
                foo.close()
            #   Now asking for images to align with the Reference image
                imgINP=[]
                foo2=open('shifts.in','w')
                imgINP=map(int, raw_input('Enter the Sl.No#s of all Images to shift and align with Ref image[ Eg: 1,2,3 ] : ').split(','))
                i=len(imgINP)
                for j in imgINP :
                    if j >= len(images) : 
                        print('Wrong serial number given.. exiting the align task..')
                        break
                    iraf.display(images[j],1)
                    print ('Press _a_ over the first star you selected in Reference image \n')
                    imx=iraf.imexam(Stdout=1)
                    Xin=eval(imx[2].split()[0])    #Using the first star as ref star for corse shift calculation
                    Yin=eval(imx[2].split()[1])
                    foo2.write(str(Xref-Xin)+'   '+str(Yref-Yin)+'\n')
                foo2.close()

                if i < 1 : break    # No images to align
                inVAR=""#Refimage
                outVAR=""# "s"+Refimage
                imgs2align=open('imgs2align.list','w')
                imgs2alignOUT=open('imgs2alignOUT.list','w')
                for j in imgINP:
                    imgs2align.write(images[j]+'\n')
                    imgs2alignOUT.write('s'+images[j]+'\n')
 #                   inVAR=inVAR+','+images[j]
 #                   outVAR=outVAR+',s'+images[j]
                imgs2align.close()
                imgs2alignOUT.close()
#                inVAR=inVAR[1:]    # To remove leading comma
#                outVAR=outVAR[1:]
                inVAR="@imgs2align.list"
                outVAR="@imgs2alignOUT.list"

                try :
                    iraf.imalign(input=inVAR, reference=Refimage, coords=Refimage+".coo", output=outVAR, shifts="shifts.in", interp_type="nearest",boundary_type="constant",trimimages="no")
                    ask=raw_input('Do you want to combine these images? [Enter _y_ for yes]: ')
                    if ask == "y" : 
                        name=raw_input('Enter name for the combined image without .fits extension: ')
                        iraf.imcombine(input=Refimage+','+outVAR, output=name+'.fits',combine="average")
                    else : print ("No combining done..\n")
                except iraf.IrafError, e :
                    print ('IRAF ERROR : Some image might be having problem. Remove it and try')
                    print e
                    print('-'*60)
                    traceback.print_exc(file=sys.stdout)
                    print('-'*60)

            
            
        print ('Moving to next filter directory\n')
        iraf.cd('../')
    print ("Moving to Next Night-------------\n")
    iraf.cd('../')


dirs.close()
print ("All alligning and combining of Images Done \n Now move required files to final directory \n And start photometry script \n Enjoy !!!----------- indiajoe@gmail.com \n")
