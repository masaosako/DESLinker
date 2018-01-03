import sys
sys.path.insert(0, './linkerCode/')

import numpy as np
import pandas as pd
import argparse
from LinkerLib import Detection
from LinkerLib import printPercentage
import time
import pickle

#takes a list of detections and returns all detections in an exposure
def findDetections(exposure, detList):
    detRes = []
    for det in detList:
        if(det.expnum == exposure):
            detRes.append(det)
    return detRes

#finds a detection with a predicted date within an exposure and attempt to grow the triplet
def growTrip(triplet, detPred, exposure, csvFile):
    detlist = wrapDets(csvFile)
    potDets = findDetections(exposure, detList)
    potTrips = []
    for det in potDets:
        newTrip = Triplet(triplet.dets.append(det))
        elements, errs = newTrip.calcOrbit()
        if(elements['a'] > 2 and elements['e'] < 1):
            potTrips.append(newTrip)
    return potTrips

#takes a csvfile of detections and wraps them with Detection class
def wrapDets(csvFile):
    df = pd.read_csv(csvFile)
    df = df.rename(columns={'SNOBJID': 'OBJID', 
            'SNFAKE_ID': 'FAKEID', 'CCDNUM': 'CCD'})
    df.columns = df.columns.str.lower()
    size = len(df['objid'])
    raList = df['ra'].tolist()
    decList = df['dec'].tolist()
    mjdList = df['mjd'].tolist()
    if 'mag' in df.columns:
        fluxList = df['mag'].tolist()
        fluxList = [10 ** ((31.4 - x) / 2.5) for x in fluxList]
    else:
        fluxList = df['flux'].tolist()
    objidList = df['objid'].tolist()
    expnumList = df['expnum'].tolist()
    ccdList = df['ccd'].tolist()
    bandList = df['band'].tolist()
    fakeidList = df['fakeid'].tolist()
    detlist = []
    startT = time.time()
    lookAhead = 0
    for y in range(size):
        printPercentage(y,size, time.time()-startT)
        det = Detection(float(raList[y]), float(decList[y]), float(mjdList[y]),
                    float(fluxList[y]), int(objidList[y]),
                    int(expnumList[y]), int(ccdList[y]), bandList[y],
                    lookAhead, int(fakeidList[y]))
        detlist.append(det)
    return detlist
    
def main():
    args = argparse.ArgumentParser()
    args.add_argument('csvFiles', nargs=1, help='path to csv file')
    args = args.parse_args()
    detlist = wrapDets(args.csvFiles[0])
    outfile = args.csvFiles[0].split('+')[-1].split('/')[-1].split('.')[0]+'.pickle'
    outfile = 'detList+' + outfile
    print('\npickling to ' + outfile) 
    with open(outfile, 'wb') as f:
        pickle.dump(detlist, f)

if __name__ == '__main__':
    main()


