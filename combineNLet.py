import sys
sys.path.insert(0,'linkerCode')
import numpy as np
import pickle
import argparse
import time
from Orbit import Orbit
from LinkerLib import Triplet, writeTriplets, pickleTriplets

'''
This program partitions the season in 10x10 (adjustable)
and combines nlet to make (n+1)let
'''

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
def combine(arr, numdet):
    # we want the merged one to contain n+1 detections
    target = numdet + 1
    combined = []
    unchanged = []
    
    minRA = arr[0,0].raLo
    minDec = arr[0,0].decLo
    stepRA = arr[0,0].raHi - minRA
    stepDec = arr[0,0].decHi - minDec

    newobjid = -1 # to be updated for the extra detection to be added

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
                    #for samp in potential:
                    #    samp.sortByMjd()
                    #toadd = [tuple([x.objid for x in y.dets]) for y in potential]
                    #if len(toadd) != len(set(toadd)):
                    #    print(str(len(toadd) - len(set(toadd))) + ' extras')
                    combined = combined + potential
                for radec in coord:
                    i = int((radec[0] - minRA) / stepRA)
                    j = int((radec[1] - minDec) / stepDec)
                    # if you want to keep the nlet, uncomment the line below
                    # arr[i,j].inactive.append(nlet1)
                    arr[i,j].nlets.remove(nlet1)
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
def checkGoodOrbit(nlets, chiSqCut = 5):
    goodnlet = []
    for nlet in nlets:
        elements, err = nlet.calcOrbit()
        chisq = nlet.getChiSq()
        if elements['a'] > 2 and elements['e'] < 1 and chisq < chiSqCut:
            goodnlet.append(nlet)
    return goodnlet

def checkDuplicates(nlets, name):
    for nlet in nlets:
        nlet.sortByMjd()
    objids = [tuple([y.objid for y in x.dets]) for x in nlets]
    print('Number of ' + name + ': ' + str(len(objids)))
    print('Number of unique ' + name + ': ' +  str(len(set(objids))))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('nlets', nargs=1, help='path to nlet pickle file')
    parser.add_argument('chisq', nargs='?', default = 5,
                         help='chi-square threshold for orbit fitter')
    args = parser.parse_args()
    nlets = pickle.load(open(args.nlets[0], 'rb'))
    
    # identify the number of detections in each Triplet object
    objcount = len(nlets[0].dets)

    t0 = time.time()
    # First get bounds of season
    # Get the array of partitioned season
    # Each cell of the array holds Triplet object with detections that fall in that region
    # if Triplet object have 5 detections (quint) then it will belong to up to 5 cells
    # Finally combine function merges nlet to form (n+1)let
    minRA, maxRA, minDec, maxDec = getBounds(nlets)
    arr = partition(nlets, minRA, maxRA, minDec, maxDec, 800, 800)
    nPlusOne, unchanged = combine(arr, objcount)
    
    # remove duplicates
    nPlusOne = removeDuplicates(nPlusOne)
    unchanged = removeDuplicates(unchanged)

    # remove same exposure
    nPlusOne = removeSameExposure(nPlusOne)
    unchanged = removeSameExposure(unchanged)

    # check for good orbits
    print('Orbit fit test in progress')
    goodNPlusOne = checkGoodOrbit(nPlusOne,float(args.chisq))
    goodUnchanged = checkGoodOrbit(unchanged, float(args.chisq))

    t1 = time.time() - t0
    print('Time taken: ' + str(t1))
    
    # dictionary for savename
    numobj = {3:'triplet', 4: 'quad', 5:'quint', 6:'hex', 7:'sept', 8:'oct'}

    # Now saving (n+1)let
    names = {3:'goodtriplet', 4:'goodquad', 5:'goodquint', 6:'goodhex',
             7:'goodsept', 8:'goodoct'}
    
    # some calculations to get stats
    total = len(goodNPlusOne) + len(goodUnchanged)
    correctUnchanged = [x for x in goodUnchanged if x.sameFake()]
    correctNPlusOne = [x for x in goodNPlusOne if x.sameFake()]
    notsame = total - len(correctUnchanged) - len(correctNPlusOne)
    inputFakeid = list(set([y.fakeid for y in x.dets for x in nlets]))
    unchangedFakeid = list(set([x.dets[0].fakeid for x in correctUnchanged]))
    nPlusOneFakeid = list(set([x.dets[0].fakeid for x in correctNPlusOne]))

    '''
    A fake with fakeid X is considered recoverable if there exists an nlet composed of
    detections only of fakeid X (i.e. all detections have fakeid X)
    '''

    print('Number of input ' + numobj[objcount] + ': ' +  str(len(nlets)))
    print('Number of input fakeids ' + str(len(inputFakeid)))
    print('Number of good fits remaining: ' + str(total))
    print('Number of good remaining ' + numobj[objcount] + ': ' + str(len(correctUnchanged)))
    print('Number of good remaining ' + numobj[objcount + 1] + ': ' +  str(len(correctNPlusOne)))
    print('Number of recoverable fakes in remaining ' + numobj[objcount] + ': ' + str(len(unchangedFakeid)))
    print('Number of recoverable fakes in remaining ' + numobj[objcount+1] + ': ' + str(len(nPlusOneFakeid)))
    print('Number of recoverable fakes remaining: ' + str(len(set(unchangedFakeid + nPlusOneFakeid))))
    
    print('\n Checking duplicates:')
    print('Input:')
    checkDuplicates(nlets, numobj[objcount])
    print('Processed:')
    checkDuplicates(goodNPlusOne, numobj[objcount+1])
    checkDuplicates(goodUnchanged, numobj[objcount])
    #print('Number of duplicates in ' + numobj[objcount+1] + ' :' + str(numdup_np1))



    savename = args.nlets[0].split('+')[-1].split('.')[0]
    saveas_n1 = names[objcount + 1] + str(args.chisq) + '+' + savename
    saveas_u = names[objcount] + 'U' + str(args.chisq) + '+' + savename
    
    # save unchanged triplet
    writeTriplets(goodUnchanged, saveas_u + '.txt', False)
    pickleTriplets(goodUnchanged, saveas_u + '.pickle')

    writeTriplets(goodNPlusOne, saveas_n1 + '.txt', False)
    pickleTriplets(goodNPlusOne, saveas_n1 + '.pickle')
    

if __name__ == '__main__':
    main()

