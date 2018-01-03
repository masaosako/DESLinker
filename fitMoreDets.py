import sys
sys.path.insert(0, './linkerCode/')

from LinkerLib import expDictionary
from LinkerLib import Detection
from LinkerLib import Triplet
from LinkerLib import printPercentage
from LinkerLib import pickleTriplets, writeTriplets
from generateDetections import wrapDets

import matplotlib.patches as patches
import numpy as np
import time
import pickle
import argparse

def withinEllipse(det, predDet):
    # The ellipse
    g_ell_center = (predDet.ra, predDet.dec)
    g_ell_width = 2*predDet.erra
    g_ell_height = 2*predDet.errb
    angle = predDet.errpa
    
    g_ellipse = patches.Ellipse(g_ell_center, g_ell_width, 
                        g_ell_height, angle=angle, fill=False, 
                        edgecolor='green', linewidth=2)
    
    cos_angle = np.cos(np.radians(180.-angle))
    sin_angle = np.sin(np.radians(180.-angle))
    
    x = det.ra
    y = det.dec
    xc = x - g_ell_center[0]
    yc = y - g_ell_center[1]

    xct = xc * cos_angle - yc * sin_angle
    yct = xc * sin_angle + yc * cos_angle 
    rad_cc = (xct**2/(g_ell_width/2.)**2) + (yct**2/(g_ell_height/2.)**2)
    
    return (rad_cc <= 1.)

#checks whetehr a single detecitons should be added
def checkDetection(nlet, predDet, newdet, chiThresh):
    debugCounter = 0
    if(withinEllipse(newdet, predDet)):
        debugCounter += 1
        newTrip = Triplet(nlet.dets[:])
        newTrip.addDetection(newdet)
        chisq = newTrip.getChiSq()
        oldsize = len(nlet.dets)
        newsize = len(newTrip.dets)
        if(chisq < chiThresh and newTrip.elements['a'] > 2
           and newsize > oldsize):
            print('\nfound new detection')
            print(newdet.toStr())
            print(newTrip.toStr())
            return newTrip, debugCounter
    return nlet, debugCounter

#adds every possible detection to the nlet
def addToNlet(nlet, expDict):
    chiThresh = 100
    counter = 0
    time0 = time.time()
    size = len(expDict)
    print('number of exposures to check: ' + str(size))
    tempLet = Triplet(nlet.dets)
    #checking every exposure
    for key, value in sorted(expDict.iteritems()):
        counter += 1
        print('\nchecking exposure ' + str(key) + 
            ' with ' + str(len(value)) + ' items')
        printPercentage(counter, size, time.time()-time0)
        mjd = value[0].mjd
        coord, erra, errb, pa = tempLet.predictPos(mjd)
        det = Detection(coord[0], coord[1], mjd, 1, 0, key, 0, 0, 0, 0) 
        det.erra = erra
        det.errb = errb
        det.pa = pa
        #checks every detection in the exposure
        debugCounter = 0 
        for pdet in value:
            tempLet, debugCount = checkDetection(tempLet, det, pdet, chiThresh)
            debugCounter += debugCount
        print('found ' + str(debugCounter) + ' within error ellipse')
    return tempLet        

def addDetections(triplet, predList, expDict):
    chiThresh = 100
    counter = 0
    time0 = time.time()
    print('number of predictions to check: ' + str(len(predList)))
    for det in predList:
        counter += 1
        printPercentage(counter, len(predList), time.time()-time0)
        potDets = []
        try:
            potDets = expDict[det.expnum]
        except KeyError:
            pass
        count2 = 0
        size2 = len(potDets)
        time1 = time.time()
        #print()
        for pdet in potDets:
            #printPercentage(count2, size2, time.time()-time1)
            count2 += 1 
            #print(str(count2) + ' of ' + str(size2) + ' checked.')
            #print(pdet.toStr())
            triplet = checkDetection(triplet, det, pdet)
            '''
            if(withinEllipse(pdet, det)):
                #print('in')
                pdets = triplet.dets[:]
                newTrip = Triplet(triplet.dets[:])
                newTrip.addDetection(pdet)
                chisq = newTrip.getChiSq()
                oldsize = len(triplet.dets)
                newsize = len(newTrip.dets)
                #print(newTrip.toStr())
                #if(pdet.fakeid == pdets[0].fakeid):
                    #time.sleep(1)
                if(chisq < chiThresh and newTrip.elements['a'] > 2 
                        and newsize > oldsize):
                    print('\nfound new detection')
                    print(pdet.toStr())
                    triplet = newTrip
                    print(triplet.toStr())
            '''
    return triplet

def main():
    args = argparse.ArgumentParser()
    args.add_argument('detList', nargs=1, help='path to csv list of detections')
    args.add_argument('predList', nargs=1, help='path to list of predicted detections')
    args = args.parse_args()

    saveName = args.predList[0].split('+')[-1].split('.')[0]
#    print(args.detList[0].split('+'))
    saveName2 = args.predList[0].split('+')[-2]
    print(saveName2) 
    print('wrapping detections: ' + args.detList[0])
    detList = wrapDets(args.detList[0])
    print('making dictionary')
    expDict = expDictionary(detList)

    print('loading prediction list: ' + args.predList[0])
    predList = pickle.load(open(args.predList[0]))
    newTripList = []
    size = len(predList)
    counter = 0
    for tup in predList:
        counter += 1
        print('*****adding to Nlet #' + str(counter) + ' of ' + str(size) + ' *****')
        print('OLD TRIPLET:')
        print(tup[0].toStr())
        
        newTriplet = addDetections(tup[0], tup[1], expDict)
        print(newTriplet.toStr())
        newTripList.append(newTriplet)
    
    name = 'crossCampaignTriplets+' + saveName2 + '+' + saveName
    writeTriplets(newTripList, name + '.txt', True)
    pickleTriplets(newTripList, name + '.pickle') 
    

if __name__=='__main__':
    main()
