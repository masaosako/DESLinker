import sys
sys.path.insert(0, 'linkerCode/')
import pickle
import argparse
import time
from combineNLet import getBounds, insertNLet, partition, combine, removeDuplicates, removeSameExposure, checkGoodOrbit, getSameRegNLet
from LinkerLib import Triplet, writeTriplets, pickleTriplets

'''
1. getBounds(nlets)
2. partition(nlets, minRA, maxRA, minDec, maxDec, 700, 700) to get array
3. combine(arr, numdet)
4. removeDuplicates(combined)
5. removeSameExposure(combined)
'''

# takes in an array of nlets and a threshold number of objects in common for it to be combined
# returns two lists:
# 1) a list of nlets that changed (got combined with some other nlet)
# 2) a list of nlets that didn't change
def generalCombine(arr, similarity, chisq):
    combined = []
    unchanged = []

    minRA = arr[0,0].raLo
    minDec = arr[0,0].decLo
    stepRA = arr[0,0].raHi - minRA
    stepDec = arr[0,0].decHi - minDec

    for x in arr:
        for y in x:
            # y is one of TripRegion
            for nlet1 in y.nlets:
                potential = []
                objid = {}
                for det1 in nlet1.dets:
                    objid[det1.objid] = det1
                # doesn't assume that all nlets have the same number of detections
                # if all nlets do have the same number of detections this could be
                # moved outside of for-loop
                numobjid = len(objid)
                othernlets, coord = getSameRegNLet(nlet1, arr, minRA, minDec, stepRA, stepDec)
                for nlet2 in othernlets:
                    newdict = objid.copy()
                    for det2 in nlet2.dets:
                        newdict[det2.objid] = det2
                    numadded = len(newdict) - numobjid
                    if numadded >= similarity:
                        potential.append(Triplet(newdict.values()))
                potential = checkGoodOrbit(potential, chisq)
                if len(potential) > 0:
                    combined = combined + potential
                else:
                    unchanged.append(nlet1)
                for radec in coord:
                    i = int((radec[0] - minRA) / stepRA)
                    j = int((radec[1] - minDec) / stepDec)
                    arr[i,j].inactive.append(nlet1)
                    arr[i,j].nlets.remove(nlet1)

    return combined, unchanged


def renewArray(arr):
    for x in arr:
        for y in x:
            y.nlets = y.nlets + y.inactive
            y.inactive = []

# compartmentalizes a list of nlets for varying n
# returns a list of list of nlets--sorted by number of detections
def compartmentalize(nlets):
    mixlets = {}
    for nlet in nlets:
        numdet = len(nlet.dets)
        if numdet in mixlets:
            mixlets[numdet].append(nlet)
        else:
            mixlets[numdet] = [nlet]
    nletList = mixlets.values()
    nletList.sort(key = lambda x: len(x[0].dets))
    return nletList
    


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('nlets', nargs=1, help='path to nlet pickle file')
    parser.add_argument('chisq', nargs='?', default = 5,
                            help='chisq threshold for good orbit')
    args = parser.parse_args()
    
    nlets = pickle.load(open(args.nlets[0], 'rb'))
    savename = args.nlets[0].split('+')[-1].split('.')[0]
    
    chisq = float(args.chisq)
    merged = []

    print('Recursively merging...')

    
    minRA, maxRA, minDec, maxDec = getBounds(nlets)
    
    while len(nlets) != 0:
        t0 = time.time()
        numdet = len(nlets[0].dets)
        print('Processing ' + str(len(nlets)) + ' ' + str(numdet) +'-lets')
        arr = partition(nlets, minRA, maxRA, minDec, maxDec, 700, 700)
        nPlusOne, unchanged = combine(arr, numdet)
        
        # remove duplicates
        nPlusOne = removeDuplicates(nPlusOne)
        unchanged = removeDuplicates(unchanged)

        # remove same exposure
        nPlusOne = removeSameExposure(nPlusOne)
        unchange = removeSameExposure(unchanged)

        # check for good orbits
        nPlusOne = checkGoodOrbit(nPlusOne, float(args.chisq))
        unchanged = checkGoodOrbit(unchanged, float(args.chisq))
        
        merged.append(unchanged)

        nlets = nPlusOne
        print('Time taken this cycle: ' + str(time.time() - t0))
    print('Done N-1 Merging')
    print('Starting Half Merge...')

    stepRA = arr[0][0].raHi - arr[0][0].raLo
    stepDec = arr[0][0].decHi - arr[0][0].decLo

    # start with nlet with most number of detections
    merged.reverse()
    remaining = []
    
    leftover = merged[0]
    
    for i in range(len(merged) - 1):
        if len(leftover) == 0:
            print('ERROR: No nlets left in this cycle')
            continue
        arr = partition(leftover, minRA, maxRA, minDec, maxDec, 700, 700)
        numdet = len(leftover[0].dets)
        print("Self merging " + str(len(leftover)) + ' ' + str(numdet) + "-let")
        sameCombined, sameUnchanged = generalCombine(arr, (numdet + 1) / 2, chisq)
        renewArray(arr)
        for oneless in merged[i+1]:
            insertNLet(oneless, arr, minRA, minDec, stepRA, stepDec)
        
        print("Second stage merging with " + str(numdet - 1) + "-let")
        mixCombined, mixUnchanged = generalCombine(arr, (len(merged[i+1][0].dets) + 1) / 2, chisq)
        remaining = remaining + sameCombined + sameUnchanged + mixCombined

        leftover = mixUnchanged

    print('Removing duplicates and same exposures')
    # remove duplicates and same exposure
    remaining = removeDuplicates(remaining)
    remaining = removeSameExposure(remaining)

    remaining = compartmentalize(remaining)
    savename = args.nlets[0].split('+')[-1].split('.')[0]
   
    for rem in remaining:
        saveas = 'good' + str(len(rem[0].dets)) + "letHM" + str(chisq) + '+' + savename
        writeTriplets(rem, saveas + '.txt', False)
        pickleTriplets(rem, saveas + '.pickle')
    print('Done')
    

if __name__ == '__main__':
    main()

