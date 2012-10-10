#!/usr/bin/env python
#USAGE:  ./SelectSource_ColorCorrectedMag.py qphot=1,2,3 source=1 filter=75 stds4beta=1,2,3 
import sys
Qphot=[]
Source=[]
Stds=[]
StdStarsMag=[[14.209,14.873,15.66,16.93],[13.679,15.312,16.920,0],[15.262,16.291,17.298,0],[15.577,16.879,18.269,0],[14.738,16.270,17.922,0],[14.131,15.946,17.893,0]]
FilterALL=['6','75','76','11']
Filter=list(FilterALL) # copy as default
# Initialising the 2-D array for storing the latest I R V data lines temperarily
iline_IRV=[]
iline_IRV_backup=[]  #For pushing and storing more than 1 pending 'same filter image' to be calculated
for i in range(len(FilterALL)) : iline_IRV.append([])
for i in range(len(FilterALL)) : iline_IRV_backup.append([])
#DONE_FLAG='Y'   #DONE_FLAG = Y/N meaing, IRV mags for lines in iline_* array has been done or not
if len(sys.argv) < 2 : 
    print('USAGE: ./SelectSource_ColorCorrectedMag.py qphot=1,2,3 source=1 filter=75 stds4beta=1,2,5,6 \n qphot, source and filter is to input the source and filter whose magnitude you want,\n stds4beta is to input standard stars to use for color eqn coefficent calculation.')
    exit()
for arg in sys.argv[1:]:
    if arg[0:6]=='qphot=' : Qphot=arg[6:].split(',')
    elif arg[0:7]=='source=' : Source=arg[7:].split(',')
    elif arg[0:7]=='filter=' : Filter=arg[7:].split(',')
    elif arg[0:10]=='stds4beta=' : Stds=arg[10:].split(',')
    else : print('Unknown input :'+arg)
#print('Searching for qphot=',Qphot,'source=',Source,'stds=',Stds, 'Filter',Filter)
StartStd=22
StartSource=18
StartQphot=8

def IVsolve(ilineI,ilineV,outFilters):  # outfilters are the filter outputs it should print.
    betaVI_array=[]
    betaIi_array=[]
    gstars1=ilineI[-2].split('_')
    gstars2=ilineV[-2].split('_')
    for star in Stds :
        if (star in gstars1) and (star in gstars2) :
            pos1=gstars1.index(star)
            pos2=gstars2.index(star)
            i=float(ilineI[StartStd+pos1*4 +2])-jd_beta[ilineI[4]]
            v=float(ilineV[StartStd+pos2*4 +2])-jd_beta[ilineV[4]]
            I=StdStarsMag[int(star)-1][0]
            V=StdStarsMag[int(star)-1][2]
            betaVI_array.append( (V-I)-1.009*(v-i) )   # for HCT CCD
            betaIi_array.append( (I-i)-0.036*(V-I) )
#            output.append(iline[StartStd+pos*4 +2])
#            output.append(iline[StartStd+(len(gstars)-1)*4+ (StartStd-StartSource)+ pos*4 +2]) #Only for McNeil
    if len(betaVI_array) > 0 :
        betaVI=sum(betaVI_array)/len(betaVI_array)
        betaIi=sum(betaIi_array)/len(betaIi_array)
        output1=[ilineI[4],ilineI[0],ilineI[2]]  # JD Fitsfile Filter 
        output2=[ilineV[4],ilineV[0],ilineV[2]]  # JD Fitsfile Filter 
        for star in Qphot :
            try :
                i=float(ilineI[StartQphot+(int(star)-1)*2])-jd_beta[ilineI[4]]
                v=float(ilineV[StartQphot+(int(star)-1)*2])-jd_beta[ilineV[4]]
                VminusI= 1.009*(v-i)+betaVI      # for HCT CCD
                I= i +0.036*VminusI +betaIi
                V=I + VminusI
                output1.append(str(round(I,3)))
                output2.append(str(round(V,3)))
#            output1.append(iline[StartQphot+(int(star)-1)*2])
            except ValueError :
                output1.append('INDEF')
                output2.append('INDEF')
        for star in Source :
            try :
                i=float(ilineI[StartSource+(int(star)-1)*4 +2])-jd_beta[ilineI[4]]
                v=float(ilineV[StartSource+(int(star)-1)*4 +2])-jd_beta[ilineV[4]]
                VminusI= 1.009*(v-i)+betaVI      # for HCT CCD
                I= i +0.036*VminusI +betaIi
                V=I+ VminusI
                output1.append(str(round(I,3)))
                output2.append(str(round(V,3)))
#            output.append(iline[StartSource+(int(star)-1)*4 +2])
            except ValueError :
                output1.append('INDEF')
                output2.append('INDEF')
        if 'I' in outFilters : print(' '.join(output1)+'  B_VI='+str(round(betaVI,3))+'  B_Ii='+str(round(betaIi,3)))
        if 'V' in outFilters : print(' '.join(output2)+'  B_VI='+str(round(betaVI,3))+'  B_Ii='+str(round(betaIi,3)))
            
#      Print the required outputs
    else :
        print('Error: No common standard stars between '+ilineI[4]+' '+ilineI[0]+' and '+ilineV[4]+' '+ilineV[0])

def RVsolve(ilineR,ilineV,outFilters):  # outfilters are the filter outputs it should print.
    betaVR_array=[]
    betaVv_array=[]
    gstars1=ilineR[-2].split('_')
    gstars2=ilineV[-2].split('_')
    for star in Stds :
        if (star in gstars1) and (star in gstars2) :
            try :
                pos1=gstars1.index(star)
                pos2=gstars2.index(star)
                r=float(ilineR[StartStd+pos1*4 +2])-jd_beta[ilineR[4]]
                v=float(ilineV[StartStd+pos2*4 +2])-jd_beta[ilineV[4]]
                R=StdStarsMag[int(star)-1][1]
                V=StdStarsMag[int(star)-1][2]
                betaVR_array.append( (V-R)-1.0277*(v-r) )   # for HCT CCD
                betaVv_array.append( (V-v)-0.1068*(V-R) )
            except ValueError:
                pass
    if len(betaVR_array) > 0 :
        betaVR=sum(betaVR_array)/len(betaVR_array)
        betaVv=sum(betaVv_array)/len(betaVv_array)
        output1=[ilineR[4],ilineR[0],ilineR[2]]  # JD Fitsfile Filter 
        output2=[ilineV[4],ilineV[0],ilineV[2]]  # JD Fitsfile Filter 
        for star in Qphot :
            try :
                r=float(ilineR[StartQphot+(int(star)-1)*2])-jd_beta[ilineR[4]]
                v=float(ilineV[StartQphot+(int(star)-1)*2])-jd_beta[ilineV[4]]
                VminusR= 1.0277*(v-r)+betaVR      # for HCT CCD
                V= v +0.1068*VminusR +betaVv
                R=V - VminusR
                output1.append(str(round(R,3)))
                output2.append(str(round(V,3)))
            except ValueError :
                output1.append('INDEF')
                output2.append('INDEF')
        for star in Source :
            try :
                r=float(ilineR[StartSource+(int(star)-1)*4 +2])-jd_beta[ilineR[4]]
                v=float(ilineV[StartSource+(int(star)-1)*4 +2])-jd_beta[ilineV[4]]
                VminusR= 1.0277*(v-r)+betaVR      # for HCT CCD
                V= v +0.1068*VminusR +betaVv
                R=V - VminusR
                output1.append(str(round(R,3)))
                output2.append(str(round(V,3)))
            except ValueError :
                output1.append('INDEF')
                output2.append('INDEF')
        if 'R' in outFilters : print(' '.join(output1)+'  B_VR='+str(round(betaVR,3))+'  B_Vv='+str(round(betaVv,3)))
        if 'V' in outFilters : print(' '.join(output2)+'  B_VR='+str(round(betaVR,3))+'  B_Vv='+str(round(betaVv,3)))
            
#      Print the required outputs
    else :
        print('Error: No common standard stars between '+ilineR[4]+' '+ilineR[0]+' and '+ilineV[4]+' '+ilineV[0])


def IRsolve(ilineI,ilineR,outFilters):  # outfilters are the filter outputs it should print.
    betaRI_array=[]
    betaIi_array=[]
    gstars1=ilineI[-2].split('_')
    gstars2=ilineR[-2].split('_')
    for star in Stds :
        if (star in gstars1) and (star in gstars2) :
            pos1=gstars1.index(star)
            pos2=gstars2.index(star)
            i=float(ilineI[StartStd+pos1*4 +2])-jd_beta[ilineI[4]]
            r=float(ilineR[StartStd+pos2*4 +2])-jd_beta[ilineR[4]]
            I=StdStarsMag[int(star)-1][0]
            R=StdStarsMag[int(star)-1][1]
            betaRI_array.append( (R-I)-0.968*(r-i) )   # for HCT CCD
            betaIi_array.append( (I-i)-0.0608*(R-I) )
#            output.append(iline[StartStd+pos*4 +2])
#            output.append(iline[StartStd+(len(gstars)-1)*4+ (StartStd-StartSource)+ pos*4 +2]) #Only for McNeil
    if len(betaRI_array) > 0 :
        betaRI=sum(betaRI_array)/len(betaRI_array)
        betaIi=sum(betaIi_array)/len(betaIi_array)
        output1=[ilineI[4],ilineI[0],ilineI[2]]  # JD Fitsfile Filter 
        output2=[ilineR[4],ilineR[0],ilineR[2]]  # JD Fitsfile Filter 
        for star in Qphot :
            try :
                i=float(ilineI[StartQphot+(int(star)-1)*2])-jd_beta[ilineI[4]]
                r=float(ilineR[StartQphot+(int(star)-1)*2])-jd_beta[ilineR[4]]
                RminusI= 0.968*(r-i)+betaRI      # for HCT CCD
                I= i +0.0608*RminusI +betaIi
                R=I + RminusI
                output1.append(str(round(I,3)))
                output2.append(str(round(R,3)))
#            output1.append(iline[StartQphot+(int(star)-1)*2])
            except ValueError :
                output1.append('INDEF')
                output2.append('INDEF')
        for star in Source :
            try :
                i=float(ilineI[StartSource+(int(star)-1)*4 +2])-jd_beta[ilineI[4]]
                r=float(ilineR[StartSource+(int(star)-1)*4 +2])-jd_beta[ilineR[4]]
                RminusI= 0.968*(r-i)+betaRI      # for HCT CCD
                I= i +0.0608*RminusI +betaIi
                R=I+ RminusI
                output1.append(str(round(I,3)))
                output2.append(str(round(R,3)))
#            output.append(iline[StartSource+(int(star)-1)*4 +2])
            except ValueError :
                output1.append('INDEF')
                output2.append('INDEF')
        if 'I' in outFilters : print(' '.join(output1)+'  B_RI='+str(round(betaRI,3))+'  B_Ii='+str(round(betaIi,3)))
        if 'R' in outFilters : print(' '.join(output2)+'  B_RI='+str(round(betaRI,3))+'  B_Ii='+str(round(betaIi,3)))
            
#      Print the required outputs
    else :
        print('Error: No common standard stars between '+ilineI[4]+' '+ilineI[0]+' and '+ilineR[4]+' '+ilineR[0])
    

# First open the Delta_Values.txt and make of dictionary of JD Vs delta (jd_beta). So that we can remove the delta correction before applying the color equations.
DeltaValuesIN=open('Delta_Values.txt','r')
prevJD='-99999'
jd_beta={}
deltaArray=[]

for deltaline in DeltaValuesIN.readlines():
    deltaline=deltaline.rstrip()
    deltaline=deltaline.split()
    if len(deltaline) > 3 :  #data lines only
        if deltaline[1] != prevJD :  #New image has started
            if len(deltaArray) > 0 :  #Just to skip the 0th  night
                jd_beta[prevJD]=sum(deltaArray)/len(deltaArray)  #Average of deltas of previous night
            deltaArray=[]
            prevJD=deltaline[1]
        deltaArray.append(float(deltaline[-2]))  #appending the delta to average in end

if len(deltaArray) > 0 :  # Entering the last image 
    jd_beta[prevJD]=sum(deltaArray)/len(deltaArray)  #Average of deltas of previous night
    deltaArray=[]
DeltaValuesIN.close()

# Now going forward for beta corrections
FinalMagIN=open('FinalDiffMag_Output.txt','r')
for iline in FinalMagIN.readlines():
    iline=iline.rstrip()
    iline=iline.split()
    output=[]
    if len(iline) == 2 and iline[0][0] !='#' and iline[1] == '-'*40 :  #New Directory/Night
        for i in range(len(FilterALL)) :
            for iline2 in iline_IRV_backup[i] :    #If any filter lines remains to be caclculated of last night.Simply print relative mag of all the images.
                if iline2[2] in Filter :  # if required by user input
                    output=[iline2[4],iline2[0],iline2[2]]  # JD Fitsfile Filter 
                    gstars=iline2[-2].split('_')
                    for star in Qphot : output.append(iline2[StartQphot+(int(star)-1)*2])
                    for star in Source : output.append(iline2[StartSource+(int(star)-1)*4 +2])
                    for mag in output[3:]:
                        if mag != 'INDEF' : 
                            print(' '.join(output)+' relativeMag')
                            break

        # if DONE_FLAG=='N' :   # Mags for previous night was not calculated
        #     lineMask=[]
        #     for i in iline_IRV : 
        #         if len(i) > 0 : lineMask.append(1)
        #         else : lineMask.append(0)
                
        #     if lineMask[0:3].count(1) == 1 : # Only one filter of IRV this night. so nothing to be done. Simply print relative mag of all the images.
        #         for iline2 in iline_IRV_backup[lineMask.index(1)]+[iline_IRV[lineMask.index(1)]] :
        #             if iline2[2] in Filter :  # if required by user input
        #                 output=[iline2[4],iline2[0],iline2[2]]  # JD Fitsfile Filter 
        #                 gstars=iline2[-2].split('_')
        #                 for star in Qphot : output.append(iline2[StartQphot+(int(star)-1)*2])
        #                 for star in Source : output.append(iline2[StartSource+(int(star)-1)*4 +2])
        #                 for mag in output[3:]:
        #                     if mag != 'INDEF' : 
        #                         print(' '.join(output)+' relativeMag')
        #                         break
        #     elif lineMask[0:3].count(1) == 2 : # Only two filters were availabel for previous night
        #         if lineMask[0]==1 and lineMask[2]==1 :  # I and V
        #             outFilters=[]
        #             if 'I' in Filter : outFilters.append('I')
        #             if 'V' in Filter : outFilters.append('V')
        #             IVsolve(iline_IRV[0],iline_IRV[2],outFilters)
        #         if lineMask[1]==1 and lineMask[2]==1 :  # R and V
        #             outFilters=[]
        #             if 'R' in Filter : outFilters.append('R')
        #             if 'V' in Filter : outFilters.append('V')
        #             RVsolve(iline_IRV[1],iline_IRV[2],outFilters)
        
        # Initialising the 2-D array for storing the latest I R V data lines temperarily
        iline_IRV=[]
        iline_IRV_backup=[]
        for i in range(len(FilterALL)) : iline_IRV.append([])
        for i in range(len(FilterALL)) : iline_IRV_backup.append([])
#        DONE_FLAG='N'
        
    if len(iline) > 4 and iline[0][0] != '#' : # and iline[2] in Filter :  #Data line,Not a discarded data line,# Required filter
        if iline[2]== FilterALL[0] : # I image
            iline_IRV[0]=list(iline)  #Saving I image line
            if len(iline_IRV[1]) != 0 : # Rimage exists. So calculate using R image
                outFilters=[]
                if '6' in Filter : outFilters.append('I')
                IRsolve(iline,iline_IRV[1],outFilters)
                if len(iline_IRV_backup[1]) !=0: # ... AND there exists R lines backups to calculate;
                    outFilters=[]
                    if '75' in Filter : 
                        outFilters.append('R')
                        for Rline in iline_IRV_backup[1] :
                            IRsolve(iline,Rline,outFilters)  #Solving all R line backups
                    iline_IRV_backup[1]=[]  # Deleting all R line backups
            elif len(iline_IRV[2]) != 0 : # Else if only V line  exists 
                outFilters=[]
                if '6' in Filter : outFilters.append('I')
                IVsolve(iline,iline_IRV[2],outFilters)
                if len(iline_IRV_backup[2]) !=0: # ... AND there exists V lines backups to calculate;
                    outFilters=[]
                    if '76' in Filter : 
                        outFilters.append('V')
                        for Vline in iline_IRV_backup[2] :
                            IVsolve(iline,Vline,outFilters)  #Solving all V line backups
                    iline_IRV_backup[2]=[]  # Deleting all V line backups
            else :      #Neither R nor V images exist. So simply backup I image for future 
                iline_IRV_backup[0].append(list(iline))
                
        if iline[2]== FilterALL[1] : # R image
            iline_IRV[1]=list(iline)  #Saving R image line
            if len(iline_IRV[0]) != 0 : # I image is present so we calculate with I
                outFilters=[]
                if '75' in Filter : outFilters.append('R')                    
                IRsolve(iline_IRV[0],iline,outFilters) # First solving the R line
                if len(iline_IRV_backup[0]) != 0: # ... AND there exists I lines backups to calculate; 
                    outFilters=[]
                    if '6' in Filter : 
                        outFilters.append('I')
                        for Iline in iline_IRV_backup[0] :
                            IRsolve(Iline,iline,outFilters)  #Solving all I line backups
                    iline_IRV_backup[0]=[]  # Deleting all I line backups
            elif len(iline_IRV[2]) != 0 : # Else if only V line exists 
                outFilters=[]
                if '75' in Filter : outFilters.append('R')                    
                RVsolve(iline,iline_IRV[2],outFilters) # First solving the R line
                if len(iline_IRV_backup[2]) != 0: # ... AND there exists V lines backups to calculate; 
                    outFilters=[]
                    if '76' in Filter : 
                        outFilters.append('V')
                        for Vline in iline_IRV_backup[2] :
                            RVsolve(iline,Vline,outFilters)  #Solving all V line backups
                    iline_IRV_backup[2]=[]  # Deleting all V line backups
            else : # Neither I nor V images exist. So simply back up R image
                iline_IRV_backup[1].append(list(iline))

        if iline[2]== FilterALL[2] : # V image
            iline_IRV[2]=list(iline)  #Saving V image line
            if len(iline_IRV[1]) != 0 : # if R line exists we calculate with R
                outFilters=[]
                if '76' in Filter : outFilters.append('V')                    
                RVsolve(iline_IRV[1],iline,outFilters) # First solving the V line
                if len(iline_IRV_backup[1]) != 0: # ... AND there exists R lines backups to calculate; 
                    outFilters=[]
                    if '75' in Filter : 
                        outFilters.append('R')
                        for Rline in iline_IRV_backup[1] :
                            RVsolve(Rline,iline,outFilters)  #Solving all R line backups
                    iline_IRV_backup[1]=[]  # Deleting all R line backups

            elif len(iline_IRV[0]) != 0 : # Else if only I line exists
                outFilters=[]
                if '76' in Filter : outFilters.append('V')                    
                IVsolve(iline_IRV[0],iline,outFilters) # First solving the V line
                if len(iline_IRV_backup[0]) != 0: # ... AND there exists I lines backups to calculate; 
                    outFilters=[]
                    if '6' in Filter : 
                        outFilters.append('I')
                        for Iline in iline_IRV_backup[0] :
                            IVsolve(Iline,iline,outFilters)  #Solving all I line backups
                    iline_IRV_backup[0]=[]  # Deleting all I line backups
            else : # Neither I nor R images exist. So simply back up V image
                iline_IRV_backup[2].append(list(iline))


        # if iline[2]== FilterALL[2] : # V image
        #     iline_IRV[2]=list(iline) #Saving V image line
        #     if len(iline_IRV[1]) == 0 and  len(iline_IRV[0]) == 0: # No R or I image yet.. Then simply backup V line..
        #         iline_IRV_backup[2].append(list(iline))
        #     if len(iline_IRV[1]) != 0 : # R image is present
        #         outFilters=[]
        #         if 'V' in Filter : outFilters.append('V')
        #         RVsolve(iline_IRV[1],iline,outFilters)
        #         if len(iline_IRV_backup[1]) != 0 : # ... AND there exists back up R images to calculate;
        #             outFilters=[]
        #             if 'R' in Filter : 
        #                 outFilters.append('R')
        #                 for Rline in iline_IRV_backup[1] :
        #                     RVsolve(Rline,iline,outFilters)  #Solving all R line backups
        #             iline_IRV_backup[1]=[]  # Deleting all R line backups
        #     if len(iline_IRV[0]) != 0 : # I Image is present.
        #         outFilters=[]
        #         if len(iline_IRV[1]) == 0 : # If R image was not present Use I to calculate V mag
        #             if 'V' in Filter : outFilters.append('V')
        #             IVsolve(iline_IRV[0],iline,outFilters)
        #         if len(iline_IRV_backup[0]) != 0 : # ... AND there exists back up I images to calculate;
        #             outFilters=[]
        #             if 'I' in Filter : 
        #                 outFilters.append('I')
        #                 for Iline in iline_IRV_backup[0] :
        #                     IVsolve(Iline,iline,outFilters)  #Solving all I line backups
        #             iline_IRV_backup[0]=[]  # Deleting all I line backups

                

#For the last night, making sure no simgle images are not calculated..                    
for i in range(len(FilterALL)) :
    for iline2 in iline_IRV_backup[i] :    #If any filter lines remains to be caclculated of last night.Simply print relative mag of all the images.
        if iline2[2] in Filter :  # if required by user input
            output=[iline2[4],iline2[0],iline2[2]]  # JD Fitsfile Filter 
            gstars=iline2[-2].split('_')
            for star in Qphot : output.append(iline2[StartQphot+(int(star)-1)*2])
            for star in Source : output.append(iline2[StartSource+(int(star)-1)*4 +2])
            for mag in output[3:]:
                if mag != 'INDEF' : 
                    print(' '.join(output)+' relativeMag')
                    break


                
        # iline_IRV[FilterALL.index(iline[2])]=iline

        # output=[iline[4],iline[0],iline[2]]  # JD Fitsfile Filter 
        # gstars=iline[-2].split('_')
        # for star in Qphot :
        #     output.append(iline[StartQphot+(int(star)-1)*2])
        # for star in Source :
        #     output.append(iline[StartSource+(int(star)-1)*4 +2])
        # for star in Stds :
        #     if star in gstars:
        #         pos=gstars.index(star)
        #         output.append(iline[StartStd+pos*4 +2])
        #         output.append(iline[StartStd+(len(gstars)-1)*4+ (StartStd-StartSource)+ pos*4 +2]) #Only for McNeil
        #     else :
        #         output.append('INDEF')
        #         output.append('INDEF')     #Only for McNeil
        # #Print only if atleast one Mag is not INDEF
        # for mag in output[3:]:
        #     if mag != 'INDEF' : 
        #         print(' '.join(output))
        #         break

FinalMagIN.close()
            
        
