#!/usr/bin/env python
# For HCT images only  
# This will convolve each image with a series of gaussin with sigma 0.5,1,1.5,2,2.5,3 and do photometry.
# IMP:  Keep ds9 open and Images4Photo.in (text file containg name of images) ready
# This script will do same photometry of all the images in text file "Images4Photo.in" and it's convolved versions
# You can use the AlignCombineImages.py to do the aligning job
# After that you can do photometry of all images together.
# Keep the following scripts also in the same directory before executing this python
# ./Creating_Log_File.sh
# ./Finding_StarID_Curser_File.sh
# Also make sure no previous attemts.mag.1 files are there in directory. Clean all previous files. <VERY IMP>
# Enjoy!!!------------------------------------------------------ indiajoe@gmail.com
import os
import os.path
import shutil
import glob
import pyfits
import sys, traceback 
import math
import numpy as np

import scipy.interpolate as interp
import numpy.ma
import warnings


# Other modules inserted inside the code are
# pyraf.iraf

def Photometry():
    iraf.noao(_doprint=0)     #Loading packages noao digiohot apphot daophot
    iraf.digiphot(_doprint=0)
    iraf.apphot(_doprint=0)
    iraf.daophot(_doprint=0)
    iraf.ptools(_doprint=0)  #for txdump

    iraf.images(_doprint=0) 
    iraf.immatch(_doprint=0) #Loading for xyxymatch, geomap, geotran of coords
    iraf.imfilter(_doprint=0)  #Loading packages for convolutiong of gauss

    iraf.phot.unlearn()   #Setting everything to default
    iraf.psf.unlearn()
    iraf.allstar.unlearn()
    iraf.daopars.unlearn()
    iraf.datapar.unlearn()
    iraf.fitskypar.unlearn()
    iraf.photpar.unlearn()
    iraf.findpars.unlearn()

    try:
        imgfile=open(MotherDIR+'/Images4Photo.in','r')
    except IOError,e:
        print('Cannot open Images4Photo.in file. Run Task #1 ')
        print(e)
        print('-'*60)
        traceback.print_exc(file=sys.stdout)
        print('-'*60)
        exit
    
    # Setting flag by checking wheter the size of qphotinput.txt is Zero or not.
    if os.stat(MotherDIR+"/qphotinput.txt")[6]!=0 : QPHOT_todo='Y'
    else : QPHOT_todo='N'
        
    imgNo=0
    for imgline in imgfile.readlines() :
        imgline=imgline.rstrip()
        imgline=imgline.split()
        wdir=imgline[0].split('/')[-2]
        if wdir != os.getcwd().split('/')[-1] : #If we are not already in img dir
            # pdir=os.getcwd().split('/') #Present directory position
            # for i in range(len(pdir)-1 - pdir.index(parentdir)): #Going back to parent directory
            #     iraf.cd('../')
            iraf.cd(MotherDIR)  #Going back to parent directory
            DIRtogo="/".join(imgline[0].split('/')[:-1]) #Now going to dir of img
            iraf.cd(DIRtogo)
            foo=open(MotherDIR+'/'+OUTPUTfile,'a')    #Appending into tbe log file The begining of a new directory
            foo.write(DIRtogo+' ---------------------------------------- \n')  # '-'*40  To mark begining of a DIrectory
            foo.close()

        
        img=imgline[0].split('/')[-1]
        filterr=imgline[1]
        intime=imgline[2]
        threshold=imgline[3]
        StartUT=imgline[4]
        print(wdir, img)
        TrueSigma=4      #Will be updated later for every img. Here it should be very small quantity
        yxdim=pyfits.getdata(img).shape  #Storing the (Ymax,Xmax) of image
        leftover=glob.glob(img+'.*')
        if len(leftover) > 0 :
            os.system('mkdir -p Leftover')
            for lft in leftover :
                os.system('mv '+lft+' Leftover/ ')


        #Calling Sextracter And find coordinates
        if not os.path.isfile(MotherDIR+"/sextractor.sex") :  #Incase the parameter file is not already created
            Sextractor_subrout(img=imgline[0])

        os.system("sextractor "+img+" -c "+MotherDIR+"/sextractor.sex -PARAMETERS_NAME "+MotherDIR+"/default.param -FILTER_NAME "+MotherDIR+"/default.conv")
#        os.system("awk 'NR>7{print $3,$5,$6}' test.cat | sort -nr | head -30 | awk '{print $2,$3}' > Bright30.coo")
        os.system("awk 'NR>7{if ($7 == 0){print $3,$5,$6}}' test.cat | sort -nr | cut -d' ' -f 2,3 | head -30 > Bright30.coo")
        #Runing xyxymatch and geo map and geotran to create new SourceT.coo , GoodStarsT.coo, BlankSky.coo
        Nmatch=32
        num_lines=0
        while num_lines < 6 :       # Set the number of stars it should atlest mach here....
            os.system("rm "+img+"xymatch.out 2> /dev/null")
            Nmatch=Nmatch-2
            iraf.xyxymatch.unlearn()
            iraf.xyxymatch(input="Bright30.coo",reference=MotherDIR+"/FirstImageTop30.coo",output=img+"xymatch.out", toler=3, matching="triangles",nmatch=Nmatch)
            # Making sure atleast a few stars were matched. otherwise geoxytran will exit with error.
            os.system("gawk '{if ($1 >0){print $0}}' "+img+"xymatch.out > matchedstars.txt") #Removing headers
            num_lines = sum(1 for line in open("matchedstars.txt"))  #Counting number of lines in the xymatch output. it should be more than 15 (size of header)
            print("Number of stars Matched= "+str(num_lines))
            if Nmatch < 5 : 
                print("Failed to find the coordinates for "+img)
                exit
        
        iraf.geomap(input=img+"xymatch.out", database=img+"rtran.db", xmin=1, xmax=yxdim[1], ymin=1, ymax=yxdim[0], interactive=0)
        
        iraf.geoxytran(input=MotherDIR+"/GoodStars.coo", output=img+"GoodStarsT.coo",database=img+"rtran.db",transforms=img+"xymatch.out")
        iraf.geoxytran(input=MotherDIR+"/Source.coo", output=img+"SourceT.coo",database=img+"rtran.db",transforms=img+"xymatch.out")
        iraf.geoxytran(input=MotherDIR+"/BlankSky.coo", output=img+"BlankSky.coo",database=img+"rtran.db",transforms=img+"xymatch.out")
        if QPHOT_todo=='Y' :
            iraf.geoxytran(input=MotherDIR+"/qphotinput.txt", output=img+"qphotinput.txt",database=img+"rtran.db",transforms=img+"xymatch.out")


        # Sanity check: To remove any new coordinates calculated lying outside image in *.coo        
        coofileLIST=['GoodStarsT.coo','SourceT.coo','BlankSky.coo']
        if QPHOT_todo=='Y' : coofileLIST.append('qphotinput.txt') 
        for coofile in coofileLIST :
            fooIN=open(img+coofile,'r')
            fooOUT=open(img+coofile+'TEMP','w')
            for star in fooIN.readlines():
                if float(star.split()[0]) > 1 and float(star.split()[0]) < yxdim[1] and float(star.split()[1]) > 1 and float(star.split()[1]) < yxdim[0] : fooOUT.write(star)
                else: print(star +": Outside the image field \n")
            fooIN.close()
            fooOUT.close()
            os.system('mv  '+img+coofile+'TEMP  '+img+coofile)

        #---------------------------------
        # Due to small error in calculation of star position, we need to create a more accurate GoodStars.coo and Source.coo
        # Plus we have to remove any saturated stars also.
        imx=iraf.imexam(input=img,frame=1,use_display=0,defkey='a',imagecur=img+'GoodStarsT.coo',Stdout=1)
        foo=open(img+'GoodStars.coo','w')    #Creating good stars coords files
        starlist=" "
        i=2
        while i < len(imx) :
            if (imx[i+1].split()[4] != 'INDEF') and (float(imx[i+1].split()[4]) > 10*math.sqrt(float(imx[i+1].split()[3]))) and (float(imx[i+1].split()[4])+float(imx[i+1].split()[3])) < float(DATAMAX) : # (Peak-sky) is not INDEF and (Peak-sky) > 10*sqrt(sky) and  Peak is not saturated
                foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
                starlist=starlist+str(i/2)+"_"   #Saving the string of good stars survived.
            else : print('Discarded: '+str(i/2)+' th number star not good of '+DIRtogo+' '+img)
            i=i+2
        foo.close()

        try :
            imx=iraf.imexam(input=img,frame=1,use_display=0,defkey='a',imagecur=img+'SourceT.coo',Stdout=1)
            Xprim=eval(imx[2].split()[0])  
            Yprim=eval(imx[2].split()[1])
            foo=open(img+'Source.coo','w')    #Creating text file contiaing coords of v1647
            i=2
            while i < len(imx) :
                foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
                i=i+2
            foo.close()
        except iraf.IrafError, e :
            print('Iraf Error: going forward with first estimated Source.coo')
            shutil.copy(img+'SourceT.coo',img+'Source.coo')
    
        #---------------------------------END of recalculation of coordinates------



        #Creating all the gauss convolved images
        OriginalIMG=img
        convIMGS=[img]
        if CONVOLVEIMG != 'NO' :  # If the CONVOLVEIMG variable is not set to NO
            for si in eval(CONVOLVEIMG) :
                iraf.gauss(input=img,output=img+'_'+str(si)+'.fits',sigma=si)
                convIMGS.append(img+'_'+str(si)+'.fits')       #List of all convolved images on which we have to do photometry
        # Now the loop of  doing the photometry for all the convolved) images
        for img in convIMGS :
            imx=iraf.imexam(input=img,frame=1,use_display=0,defkey='a',imagecur=OriginalIMG+'GoodStars.coo',Stdout=1)
            print imx           #DEBUGGING-------------------------------------------------------**
            #Calculating average FWHM
            sumfwhm=[]
            i=3
            j=0
            while i < len(imx) :               
#                sumfwhm=sumfwhm+ eval(imx[i].split()[10])
                sumfwhm.append(eval(imx[i].split()[10]))
                j=j+1
                i=i+2
            #Median FWHM is
            fwhm=np.median(sumfwhm)    # sumfwhm/j
            print('Setted value of FWHM =' + str(fwhm))
            #Calculating sky mean and stdev
            imx=iraf.imexam(input=img,frame=1,use_display=0,defkey='m',imagecur=OriginalIMG+'BlankSky.coo',Stdout=1)
            print imx            #DEBUGGING-------------------------------------------------------**
            #Calculating average Sigma and mean of sky
            sumsky=0
            sumsigma=0
            i=1
            j=0
            while i < len(imx) :               
                sumsky=sumsky+ eval(imx[i].split()[2])
                sumsigma=sumsigma+ eval(imx[i].split()[4])
                j=j+1
                i=i+1
            #Average mean,sigma and datamin are
            mean=sumsky/j
            sigma=sumsigma/j
            datamin= mean - 5*max(sigma,TrueSigma)  #Setting lowerlimit to sigma from going less than TrueSigma of original img
            print('Mean sky = '+str(mean))
            print('Mean sigma = '+str(sigma))
            print('Datamin = '+str(datamin))


            foo=open(MotherDIR+'/'+OUTPUTfile,'a')    #Appending into the log file to write output of photometry
        #    foo.write('\n -------------------------------------------- \n')
            foo.write(img +'  '+str(round(fwhm,2)) + '  '+filterr+ '  '+str(intime) +'  '+str(StartUT)+'  '+ str(round(mean,3)) +'  ' + str(round(sigma,3)) +'  '+str(round(datamin,3)) + '  ')
            foo.close()
            #Now starts photometry processes....
            # First is daopar, Press :w then :q to continue if everything is fine
            psfradi= 4*fwhm +1
            # iraf.daopars.setParam(matchrad=fwhm, psfrad=psfradi, fitrad=fwhm)
            # iraf.datapar.setParam(fwhmpsf=fwhm, sigma=sigma, datamin=datamin, datamax=65000 , readnoi=4.8 , epadu=1.22 ,itime=intime, ifilter=filterr)
            # iraf.fitskypar.setParam(annulus=270, dannulu=12 )

            iraf.daopars.setParam('matchrad',fwhm)
            iraf.daopars.setParam('psfrad',psfradi)
            iraf.daopars.setParam('fitrad',fwhm)

            iraf.datapar.setParam('fwhmpsf',fwhm)
            iraf.datapar.setParam('sigma',sigma)
            iraf.datapar.setParam('datamin',datamin)
            iraf.datapar.setParam('datamax',DATAMAX)
            iraf.datapar.setParam('readnoi',READNOISE)
            iraf.datapar.setParam('epadu',EPADU)
            iraf.datapar.setParam('itime',intime)
            iraf.datapar.setParam('ifilter',filterr)

            iraf.fitskypar.setParam('annulus',eval(ANNULUS))
            iraf.fitskypar.setParam('dannulu',eval(DANNULUS))
            

            aperture = eval(APPERTURE)  # 4*fwhm
#            threshold=6                     ## Set Detection threshold in Photometry.conf---------
            iraf.photpar.setParam('apertur',aperture)
            iraf.findpars.setParam('threshold',threshold)
            if OriginalIMG == img :
                TrueSigma=sigma
                iraf.daofind(image=img,output="default",verify=VER)
            else :
                shutil.copy(OriginalIMG+'.coo.1',img+'.coo.1')
            #Going forward to do phot
            iraf.phot(image=img,coords="default",output="default",verify=VER)
            if DOPSF == 'YES' :  # IF PSF photometry has to be done...
                #Creating the imcommands file by finding Star IDs
                os.system(MotherDIR+'/Finding_StarID_Curser_File.sh ' + img +' '+OriginalIMG+'GoodStars.coo' )
                print ("Doing psf, Non-interactively.. Using Coords of good star")
                iraf.psf(image=img, pstfile="", photfile="default", psfimage="default", opstfile="default", groupfil="default", icommands='icommands.in', verify=VER)
            #    print ("Doing psf, Interactively.. Use the a a ... f w q  sequence..")
            #    iraf.psf(image=img, pstfile="", photfile="default", psfimage="default", opstfile="default", groupfil="default")

                iraf.allstar(image=img, photfile="default", psfimage="default", allstarf="default", rejfile="default", subimage="default" ,verify=VER )
                print ("Psf photometry over")
                print ("--------------------------------------")

#%// ONLY for McNeil. Interpolating the nebula at v1647 in subtracted image and doing the phot
#%//            print("Doing Rbf interpolation of the nebula at the location of Source for "+DIRtogo+" "+img)
#%//            foo=open(OriginalIMG+'Source.coo','r')       #Reading in Coordinates of nebula infected Source
#%//            Xprim,Yprim=foo.readlines()[0].split()
#%//            foo.close()
#%//            Xprim=eval(Xprim)
#%//            Yprim=eval(Yprim)
#%//            Xprimi=int(Xprim)      #Integer coordinates
#%//            Yprimi=int(Yprim)
#%//            hdulist_sub=pyfits.open(img+'.sub.1.fits')     #psf subtracted image from allstar
#%//            nebula_sub=hdulist_sub[0].data[Yprimi-25:Yprimi+25,Xprimi-25:Xprimi+25] # 51x51 
#%//            MaskedNebula=mask_circle(nebula_sub,25,25,fwhm)  #Mask of radius fwhm at source
#%//            X,Y,Z=create1dvector_XYZ(MaskedNebula)
            #Doing the 2D Radial Basis Funtion linear interpolation using scipy.interpolate.Rbf
#%//            rbfi = interp.Rbf(X,Y, Z, smooth=2,epsilon=2,function='linear') 
#%//            XI, YI = np.meshgrid(range(nebula_sub.shape[1]), range(nebula_sub.shape[0]))
#%//            ZI = rbfi(YI, XI)
#%//            hdulist_sub[0].data[Yprimi-25:Yprimi+25,Xprimi-25:Xprimi+25]=ZI
#%//            hdulist_sub.writeto(img+'_Intp.sub.fits')  #clean Nebula in subtracted image
#%//            hdulist_sub.close()
#%//            hdulist=pyfits.open(img)
#%//            hdulist[0].data[Yprimi-25:Yprimi+25,Xprimi-25:Xprimi+25]=ZI
#%//            hdulist.writeto(img+'_Intp.fits')  #Source removed version of original image
#%//            hdulist.close()
#%//            iraf.imarith(operand1=img, op='-', operand2=img+'_Intp.sub.fits', result=img+'.New.fits')


            #Doing the phot again on stars in the new image created by subtracting interpolation
#%//            iraf.datapar.setParam('datamin',-5*max(sigma,TrueSigma))
#%//            iraf.phot(image=img+'.New.fits',coords=OriginalIMG+'Source.coo',output="default",verify=VER)
#%//            iraf.phot(image=img+'.New.fits',coords=OriginalIMG+'GoodStars.coo',output="default",verify=VER)
#%//            SecondPhotresults=iraf.txdump(textfiles=img+".New.fits.mag.1",fields="XCENTER,YCENTER,MAG",expr="yes",Stdout=1)
#%//            SecondPhotresults.extend(iraf.txdump(textfiles=img+".New.fits.mag.2",fields="XCENTER,YCENTER,MAG",expr="yes",Stdout=1))

            #Doing the phot again on Source and Good stars...
            iraf.datapar.setParam('datamin',-5*max(sigma,TrueSigma))
            iraf.phot(image=img,coords=OriginalIMG+'Source.coo',output="default",verify=VER)
            iraf.phot(image=img,coords=OriginalIMG+'GoodStars.coo',output="default",verify=VER)
            SecondPhotresults=iraf.txdump(textfiles=img+".mag.2",fields="XCENTER,YCENTER,MAG",expr="yes",Stdout=1)
            SecondPhotresults.extend(iraf.txdump(textfiles=img+".mag.3",fields="XCENTER,YCENTER,MAG",expr="yes",Stdout=1))



            iraf.hedit(img, "itime", intime, add=1, ver=0)
            
#%//            iraf.hedit(img+'_Intp.fits', "itime", intime, add=1, ver=0)
#%//            iraf.hedit(img+"_Intp.sub.fits", "itime", intime, add=1, ver=0)


#            iraf.display(img+".sub.1" , 2)
#%//            print (" Now proceeding to do qphot of Nebula from the psf subtracted image")
#%//            qphotImage=img+"_Intp.fits"
            #Qphot of the source subtracted interpolated nebula
#%//            iraf.qphot(image=qphotImage , coords=OriginalIMG+"Source.coo", cbox=5, annulus=271, dannulus=12, aperture=271, exposur="itime", epadu=EPADU ,interactive=0 )
            #Qphot of the Nebula +sources
#%//            iraf.qphot(image=img, coords=OriginalIMG+"Source.coo", cbox=5, annulus=271, dannulus=12, aperture=271, exposur="itime", epadu=EPADU , interactive=0 )


            #Doing qphot at all the points in qphotinput.txt with the corresponding parameters.
            if QPHOT_todo=='Y' :  #If there exist some qphot sources
                foo=open(OriginalIMG+"qphotinput.txt",'r')
                for qphotobj in foo.readlines():
                    qphotobj=qphotobj.rstrip()
                    obj=qphotobj.split()
                    foo2=open('qphotSource.Tcoo','w')
                    foo2.write(obj[0]+'  '+obj[1])
                    foo2.close()
                    iraf.qphot(image=img , coords='qphotSource.Tcoo', cbox=5, annulus=obj[3], dannulus=obj[4], aperture=obj[2], exposur="itime", epadu=EPADU ,interactive=0 )
                
                foo.close()
            os.system(MotherDIR+'/Creating_Log_File.sh '+img+' '+OriginalIMG+'GoodStars.coo'+' '+OriginalIMG+'Source.coo'+' '+MotherDIR+'/'+OUTPUTfile ) 
            # Writing the second phot results we txdumped X Y Mag into the list before closing the line
            foo=open(MotherDIR+'/'+OUTPUTfile,'a')  
            for line in SecondPhotresults :
                foo.write(' '+line)
            foo.write(' '+starlist+' \n') # Ending this image line with star's ID.
            foo.close()
           
            print ("Photometry of "+img+" over. \n Now proceeding to next image")
            #END of the photometry of convolved images set..
        imgNo=imgNo+1
        foo=open(MotherDIR+'/'+OUTPUTfile,'a')    #Appending into the log file to write output of photometry
        foo.write('-------------------------------------------- \n')  # '-'*44  To mark end of an image
        foo.close()

    #All photometry over
    imgfile.close()
    print("Great...Photometry of all "+str(imgNo)+" images are over...")
    print("Enjoy!!! ---------------------indiajoe@gmail.com ")


def is_number(s):   # A funtion to check wheter string s is a number or not.
    try:
        float(s)
        return True
    except ValueError:
        return False


# Below two subrouting are only for interpolation of McNeil Nebula part..
#%//def mask_circle(Matrix,a,b,r):
#%//    """ Returns the Matrix with circular mask at x,y of radius r ...
#%//    Usage  :
#%//          MaskedMatrix = myfuns.mask_circle(Matrix,x,y,r)
#%//          """
#%//    y,x=np.ogrid[-a:Matrix.shape[0]-a, -b:Matrix.shape[1]-b]
#%//    mask = x*x + y*y <= r*r
#%//    mask_array = np.zeros(Matrix.shape)
#%//    mask_array[mask] = 1 
#%//    return np.ma.masked_array(Matrix, mask=mask_array)
#%//
#%//def create1dvector_XYZ(Matrix):
#%//    """ Returns the X, Y and Z 1d vectors contining X,Y Coordinates and values Z
#%//    It doesnot include the masked pixels in matrix..
#%//    Usage :
#%//         X,Y,Z=myfuns.create1dvector_XYZ(Matrix)
#%//         """
#%//    X=[]
#%//    Y=[]
#%//    Z=[]
#%//    for i in range(Matrix.shape[0]):
#%//        for j in range(Matrix.shape[1]):
#%//            with warnings.catch_warnings():
#%//                warnings.simplefilter("ignore")
#%//                if not math.isnan(float(Matrix[i,j])) :
#%//                    X.append(i)
#%//                    Y.append(j)
#%//                    Z.append(Matrix[i,j])
#%//    return np.array(X),np.array(Y),np.array(Z)

#----------------------------------------------------

def Star_sky_subrout(img=None) :
    """ Opens the image and create Source.coo, GoodStars.coo, BlankSky.coo, Polygon.coo files"""
    backupPWD=os.getcwd()
    iraf.cd(MotherDIR)  #Going back to parent directory

    if img is None : # If No img is given, then using the first image in Images4Photo.in file
        try:
            imgfile=open(MotherDIR+'/Images4Photo.in','r')
        except IOError,e:
            print('Cannot open Images4Photo.in file. Run Task #1 ')
            print(e)
            print('-'*60)
            traceback.print_exc(file=sys.stdout)
            print('-'*60)
            exit

        imgline=imgfile.readlines()[0]
        imgline=imgline.rstrip()
        img=imgline.split()[0]
        imgfile.close()

    if not ( os.path.isfile("Source.coo") and os.path.isfile("GoodStars.coo") and os.path.isfile("BlankSky.coo") )  : #If the .coo files doesn't exist already
        iraf.display(img,1)
        print ('\n For taking coordinates of Source. Press _a_ over Primary Sources (v1647Ori)')
        imx=iraf.imexam(Stdout=1)
        foo=open('Source.coo','w')    #Creating text file contiaing coords of v1647
        i=2
        while i < len(imx) :               
            foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
            i=i+2
        foo.close()
        print ('\n For taking coordinates of good stars. Press _a_ over some good stars. \n Nonsaturated among them will be used for psf fitting.')
        print ('IMP: Press coordinate of Stars in standard required order')
        imx=iraf.imexam(Stdout=1)
        foo=open('GoodStars.coo','w')    #Creating good stars coords files
        i=2
        while i < len(imx) :               
            foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
            i=i+2
        foo.close()
        shutil.copy('GoodStars.coo','GoodStars.cooORIGINAL')   #Keeping BACKUP....
        shutil.copy('Source.coo','Source.cooORIGINAL')   #Keeping BACKUP....
        print ('\n For taking coordinates of good sky. Press _x_ over blank sky areas.')
        imx=iraf.imexam(Stdout=1)
        foo=open('BlankSky.coo','w')    #Creating blank sky coords files
        i=0
        while i < len(imx) :               
            foo.write(imx[i].split()[0] +'  '+imx[i].split()[1]+'\n')
            i=i+1
        foo.close()
        print('\n Use the ds9 to Mark circles centered at locations to do qphot.')
        print('Enter the center X Y and radius of apperture for qphot annulus and dannulus for sky')
        print('Enter values space seperated in the format below. Enter "q" to exit')
        print('X   Y   Apperture  Annulus  Dannulus ')
        foo=open('qphotinput.txt','w')    #Creating the qphot partameter file
        qphot_inp="junk"
        while (qphot_inp != "q") :
            qphot_inp=raw_input("|> ")
            boolvar=True
            for i in qphot_inp.split() : boolvar = boolvar and is_number(i)
            if boolvar and (len(qphot_inp.split()) == 5) : foo.write(qphot_inp+' \n')
            elif (qphot_inp != "q") : print("Wrong Entry. Please enter properly the 5 values. q is to stop.")
        foo.close()
    #Finished all first images data collection. Now going forward..

    print("\n All required human input of coordinates taken..")
    iraf.cd(backupPWD)



def Sextractor_subrout(img=None,N=30):
    """ Calls the Sextractor and create the sextractor parameter files if it doesn't already exists. And also create coord file of the brightest N=30 number of stars."""
    N=str(N)
    backupPWD=os.getcwd()
    iraf.cd(MotherDIR)  #Going back to parent directory
    if not os.path.isfile("sextractor.sex") : #If a config file doesn't exist already
        os.system("sextractor -d > sextractor.sex")
        if os.path.isfile("/usr/share/sextractor/default.param") :
            os.system("cp /usr/share/sextractor/default.param /usr/share/sextractor/default.conv .")
        else : print ("Error: cannot find default.param (.conv) in /usr/share/sextractor/ . \n You might have installed sextracter somewhere else")
        print("Sextractor Config file sextractor.sex and default.parm and default.conv created. \n If required u can edit it before calling Photometry")

    if img is None : # If No img is given, then using the first image in Images4Photo.in file
        try:
            imgfile=open(MotherDIR+'/Images4Photo.in','r')
        except IOError,e:
            print('Cannot open Images4Photo.in file. Run Task #1 ')
            print(e)
            print('-'*60)
            traceback.print_exc(file=sys.stdout)
            print('-'*60)
            exit

        imgline=imgfile.readlines()[0]
        imgline=imgline.rstrip()
        img=imgline.split()[0]
        imgfile.close()

    os.system("sextractor "+img+" -c sextractor.sex")
    os.system("awk 'NR>7{if ($7 == 0){print $3,$5,$6}}' test.cat | sort -nr | cut -d' ' -f 2,3 | head -"+N+" > FirstImageTop"+N+".coo")
#    os.system("awk 'NR>7{print $3,$5,$6}' test.cat | sort -nr | head -"+N+" | awk '{print $2,$3}' > FirstImageTop"+N+".coo")
    print("Brightest "+N+" stars coordinates of first image created in FirstImageTop"+N+".coo")
    iraf.cd(backupPWD)

def Cosmicrays_subrout(fileinp=None,fluxratio=6):
    """ Calls cosmicrays task in imred.crutil. Input can be single images with extension .fits or text files with first column the name of images. IMP: original will be overwritten. The default value of fluxratio=6 """
    print("Warning: Cosmic Ray removed images as overwritten on the original. I Hope you have backup...")
    print("fluxratio= "+str(fluxratio))
    iraf.imred(_doprint=0)
    iraf.crutil(_doprint=0)
    if not (fileinp is None) and fileinp[-5:]==".fits" :
        imglist=[fileinp]
    elif not (fileinp is None) :
        foo=open(fileinp,'r')
        imglist=foo.readlines()
        foo.close()
    elif fileinp is None :
        try:
            foo=open(MotherDIR+'/Images4Photo.in','r')
            imglist=foo.readlines()
            foo.close()
        except IOError,e:
            print('Cannot open Images4Photo.in file. Run Task #1 ')
            print(e)
            print('-'*60)
            traceback.print_exc(file=sys.stdout)
            print('-'*60)
            exit
    # Now doing cosmicray removal for all the images
    for img in imglist :
        img=img.rstrip()
        img=img.split()[0]
        try :
            iraf.cosmicrays(input=img, output=img,fluxratio=fluxratio,interactive='no')
        except iraf.IrafError, e :
            print('Iraf Error: Cosmic removal not done on image '+img)
    
    print("Cosmic ray removal of all images over.. \n")


def Createlist_subrout():
    """ Creates the Images4Photo.in containing the image name , filter, exposure time, threshold """
    # First creating a file with just the names of images to do photometry.
    os.system("find . -maxdepth 2 -iname '*.fits' > ImageNames.txt")
    os.system("sort ImageNames.txt > ImageNames.txtT")   #Sorting the image names
    os.system("mv ImageNames.txtT ImageNames.txt")
    # Now opening each image header and creating Images4Photo.in with filter, exposure time, threshold , UT
    foo=open('ImageNames.txt','r')
    fooOUT=open('Images4Photo.in','w')
    for img in foo.readlines():
        img=img.rstrip()
        hdulist=pyfits.open(img)
        Exptime=hdulist[0].header.get(EXPTIMEHDR)
        FilterID=hdulist[0].header.get(FILTERHDR)
        UTtime=hdulist[0].header.get(UTHDR)
        hdulist.close()
        fooOUT.write(img+' '+str(FilterID)+' '+str(Exptime)+' '+str(threshold)+' '+str(UTtime)+'  \n')
    
    fooOUT.close()
    foo.close()
    print("Images4Photo.in file created. Please edit it as required. Tip: Gawk and sed may come to use")
    print("Eg: gawk '{ if ($2==76){$4=$4-1} ; if ($2==6){$4=$4+1}; if ($2==11){$4=$4-1}; print $0}' Images4Photo.in > Images4Photo.inT ; mv Images4Photo.inT Images4Photo.in ")
    #Use FILTID for fileterID 75 =R ; 76 = V ; 6 for I for HCT
    #For Mcneil reduce the V threshold to 4 using gawk later.

def FlatCorrection_subrout(FlatStatSection):
    """ Will first combine the flats by normalising (scale=mode) and then divide images with normalised flat """
    iraf.images(_doprint=0) 
    iraf.immatch(_doprint=0) 
    iraf.imutil(_doprint=0) 
    
    iraf.imcombine.unlearn()
    
    try :
        directories=open(MotherDIR+'/directories','r')
    except IOError :  #Creating a text file containg the directories to visit if it doesn't already exist
        os.system('find . -type d -maxdepth 1 -mindepth 1 | sort > directories ')
        directories=open(MotherDIR+'/directories','r')
    for direc in directories.readlines():
        direc=direc.rstrip()
        iraf.cd(MotherDIR+'/'+direc)
        for imgDir in glob.glob('image*'):
            if len(glob.glob(imgDir+'/zs*.fits')) > 0 : #images exists
                flatdir=imgDir[5:]+'flat'
                if len(glob.glob(flatdir+'/zs*.fits')) > 0 : # atleast more than 1 flat exists...
                    os.system('ls '+flatdir+'/zs*.fits > '+flatdir+'/flatimgs')
                    iraf.imcombine (input='@'+flatdir+'/flatimgs', output= flatdir+'/'+flatdir+".fits", combine="median", scale="mode", statsec=FlatStatSection)
                    statout=iraf.imstatistics(flatdir+'/'+flatdir+".fits"+FlatStatSection,fields='mode',Stdout=1)
                    mode=float(statout[1])
                    iraf.imarith(operand1=flatdir+'/'+flatdir+".fits",op="/",operand2=mode,result=flatdir+'/'+"N"+flatdir+".fits")
                    # Now dividing all the images with normalised flat..
                    os.system('ls '+imgDir+'/zs*.fits > '+imgDir+'/imglist')
                    os.system("sed 's:/:/n:g' "+imgDir+"/imglist > "+imgDir+"/Nimglist") #prefixing n to filenames
                    iraf.imarith(operand1='@'+imgDir+'/imglist',op="/",operand2=flatdir+'/'+"N"+flatdir+".fits",result='@'+imgDir+'/Nimglist')
                else : print("ALERT: No flats in "+direc+' '+flatdir+" for flat correction ")
        print("Flat correction fo "+direc+" is over")
    print("Flat fielding of all nights are over...")
    iraf.cd(MotherDIR)                    
                
def BiasSubtract_subrout():
    """ Combine the biases in Bias directory of each night and subtract it from other images """
    iraf.noao(_doprint=0)     #Loading packages 
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)

    iraf.zerocombine.unlearn()   #Setting everything to default
    
    try :
        directories=open(MotherDIR+'/directories','r')
    except IOError :  #Creating a text file containg the directories to visit if it doesn't already exist
        os.system('find . -type d -maxdepth 1 -mindepth 1 | sort > directories ')
        directories=open(MotherDIR+'/directories','r')
    for direc in directories.readlines():
        direc=direc.rstrip()
        iraf.cd(MotherDIR+'/'+direc)        
        iraf.cd("Bias/")
        os.system("ls *.fits > zeroimgs")
        iraf.zerocombine(input= "@zeroimgs", output="Zero.fits", combine="median", ccdtype="")

        for i in ["Vflat","Rflat","Iflat","Bflat","Uflat","Haflat","imageU","imageV","imageR","imageI","imageB","imageHa"]:
            iraf.cd(MotherDIR+'/'+direc)
            iraf.cd(i)
            os.system("ls *.fits > imgs")
            iraf.imarith(operand1="@imgs",op="-",operand2="../Bias/Zero.fits",result="zs//@imgs")
        print(direc+" Night's bias subtraciton over")
        # Checking which all frames flats are missing
        for i in glob.glob('image*'):
            if len(glob.glob(i+'/*.fits')) > 0 and len(glob.glob(i[5:]+'flat/*.fits')) < 1 : # images exist in this filter but no flat
                print ('Warning: No flat for images in '+direc+' '+i)
    print("Bias subtraction of all nights over...")
    iraf.cd(MotherDIR)

def Classify_Manually_subrout():
    """ This will display one image after other, and based on user input classify images to directroies """
    try :
        directories=open(MotherDIR+'/directories','r')
    except IOError :  #Creating a text file containg the directories to visit if it doesn't already exist
        os.system('find . -type d -maxdepth 1 -mindepth 1 | sort > directories ')
        directories=open(MotherDIR+'/directories','r')
    for direc in directories.readlines():
        direc=direc.rstrip()
        iraf.cd(MotherDIR+'/'+direc)        
        for i in ["Rejects","Bias","Vflat","Rflat","Iflat","Bflat","Uflat","Haflat","imageU","imageV","imageR","imageI","imageB","imageHa"] : os.makedirs(i)
        imglist=open("rawimgs",'r')
        for img in imglist.readlines():
            img=img.rstrip()
            iraf.display(img,1)
            hdulist=pyfits.open(img)
            Comment=hdulist[0].header.get(COMMENTHDR)
            Object=hdulist[0].header.get(OBJECTHDR)
            Exptime=hdulist[0].header.get(EXPTIMEHDR)
            FilterID=hdulist[0].header.get(FILTERHDR)
            
            hdulist.close()
            print(direc,img,Object,Comment,FilterID,Exptime)
            print (" Please type d  to  reject ")
            print ("             z  for biasframe")
            print ("             fv for V flat ")
            print ("             fb for B flat ")
            print ("             fr for R flat ")
            print ("             fi for I flat ")
            print ("             fu for U flat ")
            print ("             fh for Ha flat")	
            print ("             v  for V image")
            print ("             b  for B image")
            print ("             r  for R image")
            print ("             i  for I image")
            print ("             u  for U image")
            print ("             h  for Ha image")
            verdict=""
            verdict=raw_input('Enter your Option :')
            if verdict == 'd' :  shutil.move(img,"Rejects/")
            elif verdict == 'z' : shutil.move(img,"Bias/")
            elif verdict == 'fv' : shutil.move(img,"Vflat/")
            elif verdict == 'fb' : shutil.move(img,"Bflat/")
            elif verdict == 'fr' : shutil.move(img,"Rflat/")
            elif verdict == 'fi' : shutil.move(img,"Iflat/")
            elif verdict == 'fu' : shutil.move(img,"Uflat/")
            elif verdict == 'fh' : shutil.move(img,"Haflat/")
            elif verdict == 'v' : shutil.move(img,"imageV/")
            elif verdict == 'b' : shutil.move(img,"imageB/")
            elif verdict == 'r' : shutil.move(img,"imageR/")
            elif verdict == 'i' : shutil.move(img,"imageI/")
            elif verdict == 'u' : shutil.move(img,"imageU/")
            elif verdict == 'h' : shutil.move(img,"imageHa/")
            else : print(img+" : By default Moved to Rejects")
        print("Moving all the other images still in the night directory to Rejects.")
        os.system("find . -iname \*.fits -type f -maxdepth 1 -exec mv {} Rejects/ \;")
    print("Classificatioon of images are all over..")
    iraf.cd(MotherDIR)

def OverscanTrimming_subrout():
    """ Does Overscan subtraction and Trimming of photometry images (identified by size >7Mb) And it will delete raw images """
    iraf.noao(_doprint=0)     #Loading packages noao digiohot apphot daophot
    iraf.imred(_doprint=0)
    iraf.bias(_doprint=0)

    iraf.colbias.unlearn()   #Setting everything to default
    iraf.colbias.setParam('bias',COL_BIAS)
    iraf.colbias.setParam('trim',COL_TRIM)

    try :
        directories=open(MotherDIR+'/directories','r')
    except IOError :
        #Creating a text file containg the directories to visit if it doesn't already exist
        os.system('find . -type d -maxdepth 1 -mindepth 1 | sort > directories ')
        directories=open(MotherDIR+'/directories','r')
    for direc in directories.readlines():
        direc=direc.rstrip()
        iraf.cd(MotherDIR+'/'+direc)        
        os.system("find . -maxdepth 1 -size +7M -name \*.fits  -printf '%P\n' | sort > rawimgs")
        iraf.colbias(input='@rawimgs',output='t//@rawimgs',interactive='no')
        os.system("xargs -a rawimgs rm ")  #deleting Raw images
        os.system("sed -i -n 's/^/t/p' rawimgs")  # prefixing "t" to list
        print("Overscan subtraction and Trimming of "+direc+" over.")
    print(" Overscan subtraction and Trimming of all photometry images are over...")
    iraf.cd(MotherDIR)

def Backup_subrout():
    """ Copies all the files in present directory to the ../BACKUPDIR """
    os.system('mkdir  ../'+BACKUPDIR)
    print("Copying files to ../"+BACKUPDIR)
    os.system('cp -r * ../'+BACKUPDIR)

#-------------------------------------------------------------------------------
#-----Main Program Begins here........
#def main():
try : 
    configfile=open('Photometry.conf','r')
except IOError :
    print ("Error: Copy the Photometry.conf into this directory contianing folders of each night data, before running the script.")
    exit(1)
for con in configfile.readlines():
    con=con.rstrip()
    if len(con.split()) >= 2 :
        if con.split()[0] == "VERBOSE=" :
            VER=con.split()[1]
        elif con.split()[0] == "COL_BIAS=" :
            COL_BIAS=con.split()[1]
        elif con.split()[0] == "COL_TRIM=" :
            COL_TRIM=con.split()[1]
       
        elif con.split()[0] == "THRESHOLD=" :
            threshold=con.split()[1]
        elif con.split()[0] == "EPADU=" :
            EPADU=con.split()[1]
        elif con.split()[0] == "READNOISE=" :
            READNOISE=con.split()[1]
        elif con.split()[0] == "DATAMAX=" :
            DATAMAX=con.split()[1]

        elif con.split()[0] == "APPERTURE=" :
            APPERTURE=con.split()[1]
        elif con.split()[0] == "ANNULUS=" :
            ANNULUS=con.split()[1]
        elif con.split()[0] == "DANNULUS=" :
            DANNULUS=con.split()[1]

        elif con.split()[0] == "EXPTIME=" :
            EXPTIMEHDR=con.split()[1]
        elif con.split()[0] == "FILTER=" :
            FILTERHDR=con.split()[1]
        elif con.split()[0] == "UT=" :
            UTHDR=con.split()[1]
        elif con.split()[0] == "OBJECT=" :
            OBJECTHDR=con.split()[1]
        elif con.split()[0] == "COMMENT=" :
            COMMENTHDR=con.split()[1]

        elif con.split()[0] == "OUTPUT=" :
            OUTPUTfile=con.split()[1]
        elif con.split()[0] == "BACKUP=" :
            BACKUPDIR=con.split()[1]

        elif con.split()[0] == "CONVOLVEIMG=" :
            CONVOLVEIMG=con.split()[1]
        elif con.split()[0] == "DOPSF=" :
            DOPSF=con.split()[1]

configfile.close()
MotherDIR=os.getcwd()
#    OUTPUTfilePATH=MotherDIR+'/'+OUTPUTfile
parentdir=MotherDIR.split('/')[-1]
print("Very Very Important: Backup your files first. Don't proceed without backup.\n")
print(" ---------------- The Photometry Script --------------- \n")
print("Enter the Serial numbers (space seperated if more than one task in succession) \n")
print("0  Backup files in current directory to ../"+BACKUPDIR+"\n")
print("1  Run colbias for Overscan subtraction and Trimming of Imaging data (>7Mb) & Delete raw imgs \n")
print("2  Manually classify/reject the images by displaying one by one \n")
print("3  Bias subtraction \n")
print("4  Flat Correction \n")
print("5  Make the list of images, Images4Photo.in to do Photometry \n")
print("6  Select Stars and Sky region of the field on first image \n")
print("7  Remove Cosmic Rays on all the images in Images4Photo.in. IMP:It will OVERWRITE original images.\n")
print("8  Create Sextracter config file & coordinate output of first image in this directory \n")
print("9  Do Photometry \n")
print("--------------------------------------------------------------- \n")
todo=raw_input('Enter the list : ')
todo=todo.split()
if ("1" in todo) or ("2" in todo) or ("3" in todo) or ("4" in todo) or ("6" in todo) or ("7" in todo) or ("8" in todo) or ("9" in todo) :
    from pyraf import iraf
for task in todo :
    if task == "0" :
        Backup_subrout()
    elif task == "1" :
        OverscanTrimming_subrout()
    elif task == "2" :
        Classify_Manually_subrout()
    elif task == "3" :
        BiasSubtract_subrout()
    elif task == "4" :
        FlatCorrection_subrout('[60:1780,60:1950]')

    elif task == "5" :
        Createlist_subrout()
    elif task == "6" :
        Star_sky_subrout()
    elif task == "7" :
        Cosmicrays_subrout()
    elif task == "8" :
        Sextractor_subrout()
    elif task == "9" : 
        Photometry()
print("All tasks over....Enjoy!!!_________indiajoe@gmail.com")
            

#if __name__ == "__main__":
#    main()
