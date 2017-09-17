import sys
sys.path.insert(0, '../csvDetectionfiles')

import math
import numpy as np
import pandas as pd
import argparse
from LinkerLib import Region
from LinkerLib import Detection
from LinkerLib import printPercentage
from LinkThresh import withinMag
import pickle
import time
    
#checks if two rectangles overlap given min/max ra's and dec's
'''
input: top left corner and bottom right corner of two different rectangles
output: boolean true if overlap, false otherwise
'''
def doOverlap(minRa1, maxRa1, minDec1, maxDec1,
                minRa2, maxRa2, minDec2, maxDec2):
    if(minRa1 > maxRa2 or minRa2 > maxRa1):
        return False
    if(minDec1 > maxDec2 or minDec2 > maxDec1):
        return False
    return True

#link up every detection in the cone
'''
input: a detection and a list of regions in a season
output: none, but every detection satisfying certain requirements
        is appended to the list of linked objects
'''
def linkDetection(detection, regions):
    minRa, maxRa, minDec, maxDec = detection.bounds()
    expnum = detection.expnum
    mjd = detection.mjd
    for reg in regions:
        if (doOverlap(minRa, maxRa, minDec, maxDec,
               reg.raLo, reg.raHi, reg.decLo, reg.decHi)):
            for d in reg.detections:
                if (d.expnum > expnum and d.mjd - mjd > 0.04 and
                    detection.withinCone(d)):
                    detection.linkedList.append(d)

def getSubregion(det, regions, boundary, step):
    subregion = []
    numcol = len(regions)
    numrow = len(regions[0])
    minRa, maxRa, minDec, maxDec = det.bounds()
    ((j1,j2),(i1,i2)) = getRegionIndex(minRa, minDec, maxRa, maxDec, boundary, step, regions)
    for i in range(int(i1), int(i2)):
        for j in range(int(j1), int(j2)):
                subregion.append(regions[i][j])
    return subregion

"""
Function returns indices of regions (2d array) that correspond to corners
of detection bounds
"""
def getRegionIndex(ra1, dec1, ra2, dec2, boundary, step, regions):
    stepRa = step[0]
    stepDec = step[1]
    maxRa = boundary['maxRa']
    minRa = boundary['minRa']
    maxDec = boundary['maxDec']
    minDec = boundary['minDec']

    """
    ra1 and dec1 is smaller than ra2 and dec2 respectively.
    Items from boundary are from the corners of entire region(entire 2d array).
    ra1, dec1, ra2, dec2 points are the two corners for detection bounds

    If statement for corners outside the season bounds. In this case, 
    all regions are checked. Else statement for corners within season bounds.
    Returns index of regions encompassing the corners of detection bounds
    """

    if ra1 < minRa or ra2 >= maxRa or dec1 < minDec or dec2 >= maxDec:
        return ((0, len(regions)), (0, len(regions[0])))
    else:
        delRa1 = ra1 - minRa #Find distance between ra1 and dec1
        ra1Steps = math.floor(delRa1 / stepRa) #Floor to find index - index must be integer
        delRa2 = ra2 - minRa
        ra2Steps = math.floor(delRa2 / stepRa)
        
        delDec1 = dec1 - minDec
        dec1Steps = math.floor(delDec1 / stepDec)
        delDec2 = dec2 - minDec
        dec2Steps = math.floor(delDec2 / stepDec)
        #print(ra1Steps, ra2Steps, dec1Steps, dec2Steps)
        return ((ra1Steps, ra2Steps + 1), (dec1Steps, dec2Steps + 1))


#links detections to every other detection
'''
input: list of regions with detections in each region
output: a list of detections with their linked detections
'''        
def linkDetections(regions):
    '''
    regions[0][-1] has largest RA and lowest Dec
    regions[-1][0] has lowest RA and highest Dec
    To summarize:
    for regions[i][j], increasing i would increase Dec
    and increasing j would increase RA
    '''
    startT = time.time()
    detectionLinks = []

    counter = len(regions)
    
    Lo = regions[0][0]
    Hi = regions[-1][-1]

    RaLo = Lo.raLo
    RaHi = Hi.raHi
    DecLo = Lo.decLo
    DecHi = Hi.decHi

    boundary = {'minRa':RaLo, 'maxRa':RaHi, 'minDec':DecLo, 'maxDec':DecHi}

    deltaRa =  RaHi - RaLo
    deltaDec = DecHi - DecLo


    # row corresponds to dec
    # column corresponds to ra
    numrow = len(regions)
    numcol = len(regions[0])

    step = (deltaRa/numcol, deltaDec/numrow)

    for regionlst in regions:
        for region in regionlst:
            counter2 = len(region.detections)
            startT2 = time.time()
            for det in region.detections:
                #printPercentage(len(region.detections)-counter2, len(region.detections), time.time()-startT2)
                subregion = getSubregion(det, regions, boundary, step) 
                linkDetection(det,subregion)
                detectionLinks.append(det)
                counter2 -= 1
            printPercentage(len(regions) - counter, len(regions), time.time() - startT)
            counter -= 1
    return detectionLinks

#writes the links to a file
'''
input: a detection, the file to be written
output: none, file is written
'''
def printDet(detection, txtfile):
    txtfile.write('\n' + detection.toStr() + '\n')
    txtfile.write('***********links:************** \n')
    for d in detection.linkedList:
        txtfile.write(d.toStr()+'\n')

def main():
    args = argparse.ArgumentParser()
    args.add_argument('regionsFile', nargs=1,
                        help='path to pickle file for regions in the season; ' + 
                            'filename has format regions+SNOBS_SEASON###_ML0#.pickle')
    args = args.parse_args()
    # load the regions that split the season

    print('loading region file')
    regions = pickle.load(open(args.regionsFile[0]))  
    orgFile = args.regionsFile[0].split('+')[-1].split('.')[0]
    # link up each detection with potential pairs
    print('linking detections')
    startTime = time.time()
    detectionLinks = linkDetections(np.array(regions))
    print('\n')
    print('Time taken: ' + str(time.time() - startTime))
    txtfile = 'detectionLinks+' + orgFile + '.txt'
    
    #write to a text file
    with open(txtfile, 'w+') as f:
        counter = len(detectionLinks)
        startT = time.time()
        for d in detectionLinks:
            printPercentage(len(detectionLinks) - counter, len(detectionLinks),
                time.time() - startT)
            printDet(d,f)
            counter -= 1
    #save list as a picle file
    saveName = 'detectionLinks+' + orgFile + '.pickle'
    with open(saveName, 'wb') as f:
        pickle.dump(detectionLinks, f)   

    print '\n'
if __name__ == '__main__':
    main()


