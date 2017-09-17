import sys
sys.path.insert(0, './linkerCode/')

from LinkerLib import Triplet
from LinkerLib import writeTriplets
from LinkerLib import pickleTriplets
from LinkerLib import printPercentage

import os
import time
import pickle
import argparse

def main():
    args = argparse.ArgumentParser()
    args.add_argument('folder', nargs=1, help='path to folder of chunks')
    args = args.parse_args()

    files = os.listdir(args.folder[0])
    tripList = []
    saveName = ''
    print('opening files in: ' + args.folder[0])
    time0 = time.time()
    counter = 0
    for f in files:
        counter+=1
        printPercentage(counter, len(files), time.time()-time0)
        if(f.split('+')[0] == 'goodtriplets' and f.split('.')[-1] == 'pickle'):
            print('\nopening: ' + args.folder[0] + f)
            trips = pickle.load(open(args.folder[0] + f, 'rb'))
            saveName = f.split('+')[-1].split('.')[0]
            for t in trips:
                tripList.append(t)
    
    
    writeTriplets(tripList, 'goodTriplets+' + saveName + '.txt', True)
    pickleTriplets(tripList, 'goodTriplets+' + saveName + '.pickle')
    
        


if __name__=='__main__':
    main()
