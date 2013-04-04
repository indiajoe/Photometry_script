#!/usr/bin/env python
#This python script reads through the Photmetry Output file to extract Information
import datetime
import numpy as np
import warnings
Photfile=open('Photometry_Output.txt','r')
countOUT=open('ImageCount.txt','w')
DiffMagOUT=open('DiffMag_Output.txt','w')
DeltaOUT=open('Delta_Values.txt','w')             #Deltafile
StdStars=[1,2,5,6] #,3,4
StdStarsMag=[[14.209,14.873,15.66,16.93],[13.679,15.312,16.920,0],[14.738,16.270,17.922,0],[14.131,15.946,17.893,0]]
#Complete table for McNeil Secondary Standards
#StdStars=[1,2,3,4,5,6] 
#StdStarsMag=[[14.209,14.873,15.66,16.93],[13.679,15.312,16.920,0],[15.262,16.291,17.298,0],[15.577,16.879,18.269,0],[14.738,16.270,17.922,0],[14.131,15.946,17.893,0]]

Filters=['6','75','76','11']
StartPos=16   # (16+1)th column onwards the Secondary standard stars X Y corrds start.
NoOfqphot=5  # Number of qphot sources. there are 5 columns of qphot in case of McNeil
MaxMag=30  #Maximum magnitude possible. Anything above this will be considered as coordinate in sanity check
#--------------------------------
DeltaCube=[[[] for j in range(len(Filters))] for i in range(len(StdStars))] #2d matrix of empty lists
FilterCount=[0,0,0,0]
countFlag=1
lineNo=0
JD=0
for iline in Photfile.readlines():
    lineNo=lineNo+1
    iline=iline.rstrip()
    iline=iline.split()
    DeltaArray=[]
    if (iline[0] == '-'*44) :  #End of an image
        countFlag=1 
        DiffMagOUT.write(' '.join(iline) +' \n')
        DeltaOUT.write(' '.join(iline) +' \n')      #Deltafile
    elif (iline[1] =='-'*40) :   #Begining of a Directory/Night
        if lineNo > 1 : 
            for i in FilterCount : countOUT.write(' '+str(i))  #Writing out the counts
            countOUT.write(' \n')
        countOUT.write(iline[0][2:]+' '+str(lineNo)+' ')
        for i in range(len(FilterCount)) : FilterCount[i]=0  #Reseting the counters to Zero
        countFlag=1
        DiffMagOUT.write(' '.join(iline) +' \n')
        DeltaOUT.write(' '.join(iline) +' \n')       #Deltafile
        #Calculation of Julian Date--------------
        DATE=iline[0][2:]
        if '-' in DATE : DATE=DATE.replace('-','')
        JD19900101=2447893
        if len(DATE) == 8 :
            diff = datetime.date(int(DATE[0:4]), int(DATE[4:6]), int(DATE[6:8])) - datetime.date(1990, 1, 1)
            JD=JD19900101+diff.days
        #JD calculation over--------------------
    else :
        if (countFlag == 1) :  #Incremeting Image Count counters
            try:
                i=Filters.index(iline[2])
            except ValueError:
                print('Unknown filter:'+iline[2])
                exit
            FilterCount[i]=FilterCount[i]+1 #Incrementing the image count
            countFlag=0
        #Diff MagCaculation............
        #First cheking wheter any star Mag is blank. If yes, insert INDEF there.
        sanity=0
        while sanity == 0 :
            sanity=1
            for i in range(8+NoOfqphot,len(iline)-1,3) : # Remaing X Y Mag; jumping 3 lands us on each X
                if iline[i+2] != 'INDEF' and float(iline[i+2]) > MaxMag : 
                    iline.insert(i+2,'INDEF')
                    sanity=0
                    print('Inserted an INDEF at '+str(i+2)+' in line of '+iline[0])
                    break
        #Passed sanity check..Going forward
        gstars=iline[-1].split('_')
        try:
            filt=Filters.index(iline[2])
        except ValueError:
            print('Unknown filter:'+iline[2])
            print('Skipping this image: '+iline[1])
            break
        FirstDelta=99     
        for i in range(len(StdStars)) :
            st=StdStars[i]
            if str(st) in gstars :
                pos=gstars.index(str(st))
                Mag=StdStarsMag[i]
                photmag=iline[StartPos+pos*3+2]
                if photmag != 'INDEF' : 
                    DeltaArray.append(Mag[filt]-float(photmag))
                    if i ==0 : FirstDelta= Mag[filt]-float(photmag)  #First standard star
                    elif FirstDelta != 99 :  
                        DeltaCube[i][filt].append(FirstDelta- (Mag[filt]-float(photmag)))  #The diff of deltas
        if len(DeltaArray) == 0 :
            print('No good Std star; Hence skipping this image: '+iline[0])
            continue
        #Updating the time with JD
        iline[4]=str(round(JD+(float(iline[4])/86400.0),4))

        delta=np.average(DeltaArray)
        DeltaOUT.write(iline[0]+' '+iline[4]+' '+str(DeltaArray)+' '+str(delta)+' '+str(np.std(DeltaArray))+' \n')  #Deltafile
        # Now calculation of Differential mag..
        for i in range(8,8+NoOfqphot) : # The qphot magniitudes
            if iline[i] != 'INDEF' :  iline[i]=str(round(float(iline[i])+delta,3))
        for i in range(8+NoOfqphot,len(iline)-1,3) : # Remaing X Y Mag; jumping 3 lands us on each X
            if iline[i+2] != 'INDEF' : iline[i+2]=str(round(float(iline[i+2])+delta,3)) 
        #Writing into the New output file
        DiffMagOUT.write(' '.join(iline)+' '+str(np.std(DeltaArray))+' \n')

for i in FilterCount : countOUT.write(' '+str(i))  #Writing out the last image counts
countOUT.write(' \n')
print('The header of ImageCount.txt is :: Directory LineNo #Filter1 #Filter2 ... ')
print('The DiffMag_Output.txt with differentuial magnitudes are calculated using stars : '+str(StdStars))
print('The median, mean and std of delta of each secondary standars used  for each filter are:')
for i in range(len(StdStars)):
    print('----- '+str(StdStars[i]))
    for j in range(len(Filters)):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            print(Filters[j],np.median(DeltaCube[i][j]),np.mean(DeltaCube[i][j]),np.std(DeltaCube[i][j]))

    
countOUT.close()
DiffMagOUT.close()
Photfile.close()
DeltaOUT.close()        #Deltafile
# For Linear fitting of fwhm vs magnitude and extrapolating back to particular fwhm value of mag.
FWHM=4  #The FWHM at which we should calculate the magnitude by fitting straight line.

DiffMagfile=open('DiffMag_Output.txt','r')
FinalMagOUT=open('FinalDiffMag_Output.txt','w')
#import numpy.linalg as la
countFlag=1
for iline in DiffMagfile.readlines():
    iline=iline.rstrip()
    iline=iline.split()
    if (iline[0] == '-'*44) :  #End of an image
        if ImgFlag == 1 :  #Only if the image is loaded in ilineCP
        #Doing the linear least square fit of mag vs fwhm
            for i in range(8,8+NoOfqphot) : # The qphot magniitudes
                X=[]
                Y=[]
                for j in range(len(ilineCP[i])) :
                    if ilineCP[i][j] != 'INDEF' : 
                        Y.append(float(ilineCP[i][j]))  #mag
                        X.append(float(ilineCP[1][j]))  #fwhm
                if len(Y) <= 1 :      #Atleast 2 points required for straight line fit
                    if len(Y) ==1: print('LSF not done for one point in '+str(i)+'th column of '+ilineCP[0])
                    ilineCP[i]='INDEF'    
                    continue
                A = np.vstack([X, np.ones(len(X))]).T
                #m,c = np.linalg.lstsq(A, Y)[0] # y=mx+c is the straight line 
                #residue = np.linalg.lstsq(A, Y)[1][0]
                if len(Y) >2 :([m,c],[residue])=np.linalg.lstsq(A,Y)[0:2] # y=mx+c is straight line & 'residue' is least square residue
                else:
                    m,c=np.linalg.lstsq(A, Y)[0]
                    residue=99        #Only 2 points, hence residue is meaning less

                ilineCP[i]=str(round(m*FWHM+c,3))+' '+str(round(residue,3))  #Replacing the array with string 'Mag residue'

            for i in range(8+NoOfqphot,len(ilineCP)-2,3) : # Remaing X Y Mag; jumping 3 lands us on each X
                X=[]
                Y=[]
                for j in range(len(ilineCP[i+2])) :
                    if ilineCP[i+2][j] != 'INDEF' : 
                        Y.append(float(ilineCP[i+2][j]))  #mag
                        X.append(float(ilineCP[1][j]))  #fwhm
                if len(Y) <= 1 :      #Atleast 2 points required for straight line fit
                    if len(Y) ==1: print('LSF not done for one point in '+str(i)+'th column of '+ilineCP[0])
                    ilineCP[i+2]='INDEF INDEF'    
                    continue
                A = np.vstack([X, np.ones(len(X))]).T
                if len(Y) > 2: ([m,c],[residue])=np.linalg.lstsq(A,Y)[0:2] # y=mx+c is the straight line and 'residue' is lsf residue
                else :
                    m,c=np.linalg.lstsq(A, Y)[0]
                    residue=99        #Only 2 points, hence residue is meaning less
                ilineCP[i+2]=str(round(m*FWHM+c,3))+' '+str(round(residue,3))  #Replacing the array with string 'Mag residue'

            ilineCP[1]=str(FWHM)  #Replacing array with string of fwhm
            FinalMagOUT.write(' '.join(ilineCP) +' \n')
            ImgFlag=0  # Unloaded the image from ilineCP
            countFlag=1 

    elif (iline[1] =='-'*40) :   #Begining of a Directory/Night
        countFlag=1
        FinalMagOUT.write(' '.join(iline) +' \n')
    else :
        if (countFlag == 1) :  #First line of an Image
            ilineCP=iline[:] # create a fresh copy
            ImgFlag=1    #Image loaded Flag set to 1
            #Converting the Mag colums into a python List
            ilineCP[1]=[ilineCP[1]]  #fwhm
            for i in range(8,8+NoOfqphot) : # The qphot magniitudes
                ilineCP[i]=[ilineCP[i]]
            for i in range(8+NoOfqphot,len(iline)-2,3) : # Remaing X Y Mag; jumping 3 lands us on each X
                ilineCP[i+2]=[ilineCP[i+2]]
            countFlag=0
        else :  # Gaussian Convolved images
            ilineCP[1].append(iline[1])  # appending fwhm
            for i in range(8,8+NoOfqphot) : # appending qphot magniitudes
                ilineCP[i].append(iline[i])
            for i in range(8+NoOfqphot,len(iline)-2,3) : # appending Remaing X Y Mag; jumping 3 lands us on each X
                ilineCP[i+2].append(iline[i+2])


print('FinalDiffMag_Output.txt created with Least square fitted Magnitude and residue of lsq')
FinalMagOUT.close()
DiffMagfile.close()

