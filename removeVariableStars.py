import sys
sys.path.insert(0,'./linkerCode/')
import numpy as np
import argparse
import time
import pandas as pd
from scipy.spatial import KDTree
from LinkerLib import printPercentage

#distance between two points
def dist(pt1, pt2):
    deltax = pt1[0]-pt2[0]
    deltay = pt1[1]-pt2[1]
    return (deltax*deltax + deltay*deltay)

#takes two tuples and returns whether their values are within floating points of each other
#takes the last element each tuple and returns true if they're within floating point
def checkWithinFloating(props1, props2, thresh = 0.0000000001):
    l1 = len(props1)
    l2 = len(props2)
    if(not l1 == l2):
        return False, False
    withinFloat = True
    withinSorted = True
    for x in range(l1-2):
        val1 = float(props1[x])
        val2 = float(props2[x])
        if(abs(val1-val2) > thresh):
            withinFloat = False
    if(not props1[-2] == props2[-2]):
        withinFloat = False
    if(abs(float(props1[-1])-float(props2[-1])) > thresh):
        withinSorted = False
    return withinFloat, withinSorted

#takes in a dataframe and removes duplicate detections
def rmDuplicates(detDF, sortedInd):
    print('size before removing duplicates: ' + str(len(detDF))) 
    rmInd = []
    expList = detDF['EXPNUM'].tolist()
    ccdList = detDF['CCDNUM'].tolist()
    raList = detDF['RA'].tolist()
    decList = detDF['DEC'].tolist()
    # temporary fix for mag/flux issue
    if 'MAG' in detDF.columns:
        fluxList = detDF['MAG'].tolist()
        fluxList = [10 ** ((31.4 - x) / 2.5) for x in fluxList]
    else:
        fluxList = detDF['FLUX'].tolist()
    bandList = detDF['BAND'].tolist()
    sortedInd = detDF[sortedInd].tolist()

    objidList = detDF['OBJID'].tolist()
    time0 = time.time()
    for x in range(len(expList)):
        if x in rmInd:
            continue
        printPercentage(x, len(expList), time.time()-time0)
        props1 = (expList[x], ccdList[x], raList[x], decList[x],
                    fluxList[x], bandList[x], sortedInd[x])
        objid1 = objidList[x]
        withinSorted = True 
        while withinSorted and x<len(expList)-1:
            x += 1
            props2 = (expList[x], ccdList[x], raList[x], decList[x],
                        fluxList[x], bandList[x], sortedInd[x])
            objid2 = objidList[x]
            withinFloating, withinSorted = checkWithinFloating(props1, props2)
            if(withinFloating):
                rmInd.append(x)   
    print('\nnumber removed: ' + str(len(rmInd)))
    detDF = detDF.drop(detDF.index[rmInd])
    print('size after removing: ' + str(len(detDF)))
    return detDF

#takes in dataframe of detections and dataframe of variable stars
#writes the nonvariable stars in outfile
def rmNonVariables(detDF, varDF, distThresh = 0.000027778):
    print('extracting coordinates')
    #convert columns into list of points
    varRAs = varDF['RA'].tolist()
    varDECs = varDF['DEC'].tolist()
    varPoints = zip(varRAs, varDECs)
    kdTree = KDTree(varPoints)
    detRAs = detDF['RA'].tolist()
    detDECs = detDF['DEC'].tolist()
    #remove the rows that are close to a variable star
    time0 = time.time()
    findTime = 0
    distTime = 0
    dropTime = 0
    print('removing variable stars')
    rmList = []
    print('distance threshold:' + str(distThresh))
    print('number of detections:' + str(len(detDF)))
    for x in range(len(detDF)):
        printPercentage(x, len(detDF), time.time()-time0)
        pt1 = (float(detRAs[x]), float(detDECs[x]))
        dist, i = kdTree.query(pt1, 1)
        if(dist < distThresh):
            rmList.append(x)
    print('\nnumber of detections removed: ' + str(len(rmList)))
    detDF = detDF.drop(detDF.index[rmList])
    print('number of detections remaining: ' + str(len(detDF)))
    return detDF

def main():
    args = argparse.ArgumentParser()
    args.add_argument('detections', nargs=1, help='path/filename of csv detection file; ' + 
                        'filename has format SNOBS_SEASON###_ML0#.csv')
    args.add_argument('variables', nargs=1, help='path to csv file of variable stars')
    args.add_argument('distThresh', nargs='?', default='0.00027778', 
                help='maximum distance to be considered same object (in degrees)')
    args = args.parse_args()
    savename = args.detections[0].split('+')[-1].split('/')[-1].split('.')[0]
    detDF = pd.read_csv(args.detections[0])
    varDF = pd.read_csv(args.variables[0])
    detDF = rmNonVariables(detDF, varDF, float(args.distThresh))
    detDF = rmDuplicates(detDF, 'RA')
    detDF.to_csv('varRM+' + savename + '.csv')

if __name__=='__main__':
    main()
