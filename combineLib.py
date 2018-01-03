import sys
sys.path.insert(0,'linkerCode')
import numpy as np
import pickle
import argparse
import time
from Orbit import Orbit
from LinkerLib import Triplet, writeTriplets, pickleTriplets, printPercentage

'''
This program partitions the season in 10x10 (adjustable)
and combines nlet to make (n+1)let
'''

class resizableArray:
    def __init__(self):
        self.arr = np.empty([1],dtype=object)
        self.ptr = 0

    def append(self,obj):
        if self.arr[-1] is None:
            self.arr[self.ptr] = obj
            self.ptr += 1
        else:
            #print("Updating array size to: " + str(len(self.arr) * 2))
            newArr = np.empty([len(self.arr) * 2], dtype=object)
            for i in range(len(self.arr)):
                newArr[i] = self.arr[i]
            self.arr = newArr
            self.arr[self.ptr] = obj
            self.ptr += 1

    def toList(self):
        return [x for x in list(self.arr) if x is not None]


# This object will be the cells of the array to contain the nlets
class NLetRegion:
    def __init__(self, raLo, raHi, decLo, decHi):
        self.raLo = raLo
        self.raHi = raHi
        self.decLo = decLo
        self.decHi = decHi
        self.nlets = []
        self.inactive = [] # use this if you want to keep nlets after combine()

    def add(self, nlet):
        self.nlets.append(nlet)

# returns the bounds of season
def getBounds(nlets):
    minRA = 360
    maxRA = -360
    minDec = 360
    maxDec = -360
    for nlet in nlets:
        ra = [x.ra for x in nlet.dets]
        dec = [x.dec for x in nlet.dets]
        minRA = min(minRA, min(ra))
        maxRA = max(maxRA, max(ra))
        minDec = min(minDec, min(dec))
        maxDec = max(maxDec, max(dec))
    minRA -= 0.01
    maxRA += 0.01
    minDec -= 0.01
    maxDec += 0.01
    return minRA, maxRA, minDec, maxDec

# Inserts all of the nlets into the array
# Helper function for partition()
def insertNLet(nlet, arr, minRA, minDec, stepRA, stepDec):
    coord = zip([x.ra for x in nlet.dets], [x.dec for x in nlet.dets])
    for radec in coord:
        i = int((radec[0] - minRA) / stepRA)
        j = int((radec[1] - minDec) / stepDec)
        arr[i,j].add(nlet)
    

# Partitions the season
# Returns an array with each cell containing Triplet object with detection
# that falls into that cell
def partition(nlets, minRA, maxRA, minDec, maxDec, numcol, numrow):
    # initialize array
    stepRA = (maxRA - minRA) / numcol
    stepDec = (maxDec - minDec) / numrow
    numrow = int(numrow)
    numcol = int(numcol)
    arr = np.empty((numrow, numcol), dtype = object)
    for i in range(numcol):
        for j in range(numrow):
            arr[i,j] = NLetRegion(minRA + i * stepRA, minRA + (i + 1) * stepRA,
                                    minDec + j * stepDec, minDec + (j + 1) * stepDec)
    
    # insert nlets into NLetRegion
    for nlet in nlets:
        insertNLet(nlet, arr, minRA, minDec, stepRA, stepDec)
    return arr

# returns all of the Triplet objects that belongs to common cell in the array
# and the coordinates of detections of input nlet
def getSameRegNLet(nlet, arr, minRA, minDec, stepRA, stepDec):
    coord = zip([x.ra for x in nlet.dets], [x.dec for x in nlet.dets])
    nlets = []
    for radec in coord:
        i = int((radec[0] - minRA) / stepRA)
        j = int((radec[1] - minDec) / stepDec)
        nlets = nlets + arr[i,j].nlets
    nlets = list(set(nlets))
    nlets.remove(nlet)
    return nlets, coord

# returns nlet with lowest chisq
def lowestChiSq(nlets):
    minchisq = nlets[0].getChiSq()
    lowest = nlets[0]
    for i in range(1, len(nlets)):
        chisq = nlets[i].getChiSq()
        if chisq < minchisq:
            lowest = nlets[i]
            minchisq = chisq
    return lowest
    

# takes in arr and number of detections in each Triplet object
def combine(arr, numdet, progress=False):
    # we want the merged one to contain n+1 detections
    target = numdet + 1
    #combined = []
    #unchanged = []
    
    combined = resizableArray()
    unchanged = resizableArray()

    minRA = arr[0,0].raLo
    minDec = arr[0,0].decLo
    stepRA = arr[0,0].raHi - minRA
    stepDec = arr[0,0].decHi - minDec

    newobjid = -1 # to be updated for the extra detection to be added
    if progress:
        count = 0
        time0 = time.time()
        size = 0
        for i in arr:
            for j in i:
                size += len(j.nlets)
        size /= numdet

    for x in arr:
        for y in x:
            # y is one of TripRegion
            for nlet1 in y.nlets:
                potential = []
                objid = [i.objid for i in nlet1.dets]
                addedId = []
                othernlets, coord = getSameRegNLet(nlet1, arr, minRA, minDec, stepRA, stepDec)
                for nlet2 in othernlets:
                    detgroup = [x for x in nlet1.dets]
                    for det in nlet2.dets:
                        if det.objid not in objid and det.objid not in addedId:
                            detgroup.append(det)
                            newobjid = det.objid
                    # check that there are n+1 objects
                    if len(detgroup) == target:
                        potential.append(Triplet(detgroup))
                        addedId.append(newobjid)
                        
                if len(potential) == 0:
                    unchanged.append(nlet1)
                else:
                    for pot in potential:
                        combined.append(pot)
                if progress:
                    count += 1
                    printPercentage(count, size, time.time() - time0)
                for radec in coord:
                    i = int((radec[0] - minRA) / stepRA)
                    j = int((radec[1] - minDec) / stepDec)
                    # if you want to keep the nlet, uncomment the line below
                    # arr[i,j].inactive.append(nlet1)
                    arr[i,j].nlets.remove(nlet1)
    return combined.toList(), unchanged.toList()


def generalCombine(arr, similarity, keep=False, progress=False):
    """
    arr - numpy array of NLetRegion
    similarity - float between 0 and 1
    keep - set it to True to keep nlets in inactive list of NLetRegion
    progress - set it to True for showing progress
    """
    combined = resizableArray()
    unchanged = resizableArray()

    minRA = arr[0,0].raLo
    minDec = arr[0,0].decLo
    stepRA = arr[0,0].raHi - minRA
    stepDec = arr[0,0].decHi - minDec

    if progress:
        count = 0
        time0 = time.time()
        uniqueNLets = {}
        for i in arr:
            for j in i:
                for nlet in j:
                    nlet.sortByMjd()
                    uniqueNlets[tuple([x.objid for x in nlet.dets])] = nlet
        size = len(uniqueNlets.values())

    for x in arr:
        for y in arr:
            for nlet in y:
                potential = []
                objid = set([x.objid for x in nlet.dets])
                target = int(similarity * len(nlet.dets))
                othernlets, coord = getSameRegNLet(nlet, arr, minRA, minDec,
                                                    stepRA, stepDec)
                for other in othernlets:
                    otherID = set([x.objid for x in others.dets])
                    if len(objid.intersection(otherID)) >= target:
                        detgroup = [x for x in nlet.dets]
                        potential.append(Triplet(list(set([x for x in other.dets]
                                                        + detgroup))))
                
                if len(potential) != 0:
                    for pot in potential:
                        combined.append(pot)
                else:
                    unchanged.append(nlet)

                if progress:
                    count += 1
                    printPercentage(count, size, time.time() - time0)

                for radec in coord:
                    i = int((radec[0] - minRA) / stepRA)
                    j = int((radec[1] - minDec) / stepDec)
                    if keep:
                        arr[i,j].inactive.append(nlet)
                    arr[i,j].nlets.remove(nlet)
    return combined, unchanged


# remove duplicates
def removeDuplicates(nlets):
    for nlet in nlets:
        nlet.sortByMjd()
    unique = {}
    for nlet in nlets:
        unique[tuple([x.objid for x in nlet.dets])] = nlet
    return unique.values()
        
def removeSameExposure(nlets):
    removed = []
    for nlet in nlets:
        expnum = [x.expnum for x in nlet.dets]
        if len(set(expnum)) == len(expnum):
            removed.append(nlet)
    return removed


# checks for good orbit using semi-major axis, eccentricity, chisq
def checkGoodOrbit(nlets, chiSqCut = 5, progress=False):
    goodnlet = []
    time0 = time.time()
    size = len(nlets)
    count = 0
    for nlet in nlets:
        elements, err = nlet.calcOrbit()
        chisq = nlet.getChiSq()
        if elements['a'] > 2 and elements['e'] < 1 and chisq < chiSqCut:
            goodnlet.append(nlet)
        if progress:
            count += 1
            printPercentage(count, size, time.time() - time0)
    if progress:
        print("")
    return goodnlet

def checkDuplicates(nlets, name):
    for nlet in nlets:
        nlet.sortByMjd()
    objids = [tuple([y.objid for y in x.dets]) for x in nlets]
    print('Number of ' + name + ': ' + str(len(objids)))
    print('Number of unique ' + name + ': ' +  str(len(set(objids))))

