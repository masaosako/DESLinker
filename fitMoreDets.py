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
    arcSecToDegree = 0.000277778
    g_ell_center = (predDet.ra, predDet.dec)
    g_ell_width = 2*predDet.erra*arcSecToDegree
    g_ell_height = 2*predDet.errb*arcSecToDegree
    angle = predDet.pa
    
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
    timeWrap = 0
    timeChiSq = 0
    timeA = time.time()
    inEllipse = withinEllipse(newdet, predDet)
    timeB = time.time()
    if(inEllipse):
        print('\nra:' + str(predDet.ra) + ' dec:' + str(predDet.dec))
        print('erra:' + str(predDet.erra) + ' errb:'+ str(predDet.errb) + ' pa:' + str(predDet.pa))
        print(newdet.toStr())
        debugCounter += 1
        newTrip = Triplet(nlet.dets[:])
        newTrip.addDetection(newdet)
        timeB = time.time()
        chisq = newTrip.getChiSq()
        oldsize = len(nlet.dets)
        newsize = len(newTrip.dets)
        ##
        timeC = time.time()
        timeWrap = timeB-timeA
        timeChiSq = timeC-timeB
        #print('timewrap' + str(timeWrap))
        #print('timechisq' + str(timeChiSq))
        if(chisq < chiThresh and newTrip.elements['a'] > 2
           and newsize > oldsize):
            print('\nfound new detection')
            print(newdet.toStr())
            print(newTrip.toStr())
            return newTrip, debugCounter, (timeWrap, timeChiSq)
    timeWrap = timeB - timeA
    return nlet, debugCounter, (timeWrap, timeChiSq)

def dontworryaboutwhatthisdoes(dictionaryThing):
    for x in dictionaryThing:
        print(x[0])

#adds every possible detection to the nlet
def addToNlet(nlet, expDict):
    time1 = time.time()
    expLst = sorted([x.expnum for x in nlet.dets])
    ####
    time2 = time.time()
    firstExp = expLst[len(expLst)/2]
    lastExp = firstExp
    print('first exp: ' + str(firstExp) + ' last exp: ' + str(lastExp))
    chiThresh = 100
    counter = 0
    ####
    size = len(expDict)
    print('number of exposures to check: ' + str(size))
   
    tempLet = Triplet(nlet.dets)
    
    #only look at exposures after last and before first
    time4 = time.time()
    sortedDict = sorted(expDict.iteritems())
    greater = [x for x in sortedDict if x[0]>lastExp]
    lesser = [x for x in sortedDict if x[0]<firstExp]
    lesser.reverse()
    time0 = time.time()
    #dontworryaboutwhatthisdoes(greater)
    #dontworryaboutwhatthisdoes(lesser)
    timeList = time2-time1
    timeInit = time4-time2
    timeSort = time0-time4
    timePrint = 0
    timePred = 0
    timeWrap = 0
    timeCheck = 0
    timeCheckWrap = 0
    timeCheckChi = 0
    totalCount = 0
    #checking every exposure
    for key, value in greater:
        counter += 1
        timeA = time.time()
        print '\nchecking exposure ' + str(key) + ' with ' + str(len(value)) + ' items', 
        ###    
        
        #printPercentage(counter, size, time.time()-time0)
        ###
        timeB = time.time()
        mjd = value[0].mjd
        coord, erra, errb, pa = tempLet.predictPos(mjd)
        #print('pa' + str(pa))
        ###
        timeC = time.time()
        det = Detection(coord[0], coord[1], mjd, 1, 0, key, 0, 0, 0, 0) 
        det.erra = erra
        det.errb = errb
        det.pa = pa
        print 'erra: ' + str(erra) + ' errb: ' + str(errb),
        ###
        timeD = time.time()
        #checks every detection in the exposure
        debugCounter = 0 
        for pdet in value:
            tempLet, debugCount, (timeF, timeG) = checkDetection(
                    tempLet, det, pdet, chiThresh)
            debugCounter += debugCount
            timeCheckWrap += timeF
            timeCheckChi += timeG
        timeE = time.time()
        print 'found ' + str(debugCounter) + ' within error ellipse',
        totalCount += debugCounter
        timePrint += timeB-timeA
        timePred += timeC-timeB
        timeWrap += timeD-timeC
        timeCheck += timeE-timeD
    
    for key, value in lesser:
        counter += 1
        timeA = time.time()
        print '\nchecking exposure ' + str(key) + ' with ' + str(len(value)) + ' items',
        #printPercentage(counter, size, time.time()-time0)
        ##
        timeB = time.time()
        mjd = value[0].mjd
        coord, erra, errb, pa = tempLet.predictPos(mjd)
        timeC = time.time()
        det = Detection(coord[0], coord[1], mjd, 1, 0, key, 0, 0, 0, 0) 
        det.erra = erra
        det.errb = errb
        det.pa = pa
        print 'erra: ' + str(erra) + ' errb: ' + str(errb),

        ##
        timeD = time.time()
        #checks every detection in the exposure
        debugCounter = 0 
        for pdet in value:
            tempLet, debugCount, (timeF, timeG) = checkDetection(
                tempLet, det, pdet, chiThresh)
            debugCounter += debugCount
            timeCheckWrap += timeF
            timeCheckChi += timeG
        totalCount +=debugCounter
        timeE = time.time()
        print 'found ' + str(debugCounter) + ' within error ellipse',
        timePrint += timeB-timeA
        timePred += timeC-timeB
        timeWrap += timeD-timeC
        timeCheck += timeE-timeD
    print('timeList: ' + str(timeList))
    print('timeInit: ' + str(timeInit))
    print('timeSort: ' + str(timeSort))
    print('timePrint: ' + str(timePrint))
    print('timePred: ' + str(timePred))
    print('timeWrap: ' + str(timeWrap))
    print('timeCheck: ' + str(timeCheck))
    print('timeCheckWrap: ' + str(timeCheckWrap))
    print('timeCheckChi: ' + str(timeCheckChi))
    print('total detections checked: ' + str(totalCount))
    return tempLet        

def addDetections(triplet, predList, expDict):
    chiThresh = 10
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
