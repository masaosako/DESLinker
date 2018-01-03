import sys
sys.path.insert(0, './linkerCode/')

import pickle
import pandas as pd
import argparse
import time
import numpy as np

from LinkerLib import Detection
from LinkerLib import printPercentage
from LinkerLib import Triplet
from LinkerLib import expDictionary
from LinkerLib import pickleTriplets, writeTriplets

from generateDetections import wrapDets
from detToDictionary import detPrediction
from fitMoreDets import addDetections


def main():
    args = argparse.ArgumentParser()
    #args.add_argument('csvDetList', nargs=1, help='path to csv file of list of detections')
    args.add_argument('expList', nargs=1, help='a corresponding list of exposures for that csvList')
    args.add_argument('potentialNlets', nargs=1, help='Nlets we want to add on to')
    args = args.parse_args()
    #print('\nwrapping detections ' + args.csvDetList[0])
    #detlist = wrapDets(args.csvDetList[0])
    print('\nloading nlets ' + args.potentialNlets[0])
    nlets = pickle.load(open(args.potentialNlets[0], 'rb'))
    print('\npredicting nlets in exposures')
    predDets = detPrediction(nlets, args.expList[0])
    
    newTripList = []
    #print('\nmaking dictionary')
    #expDict = expDictionary(detlist)
    print('\n dictionary finished')
    saveName = args.potentialNlets[0].split('+')[-1].split('.')[0]
    saveName2 = args.potentialNlets[0].split('+')[-2].split('/')[-1]
    size = len(predDets)
    savesize = 5
    counter = 0
    print('\nsaving predictions')
    time0 = time.time()
    while(counter < size):
        printPercentage(counter, size, time.time()-time0)
        predChunkTup = []
        for x in range(min(savesize, size - counter)):
            predChunkTup.append(predDets[counter])
            predChunkTup[x][0].orbit = 0
            counter += 1
        chunkNum = counter/savesize
        with open('listOfPredictions+Chunk' + str(chunkNum) + 
                '_' + saveName2 + '+' + saveName + '.pickle', 'wb') as f:
            pickle.dump(predChunkTup, f)   
    '''
    
    for tup in predDets:
        counter += 1
        print('*****adding to Nlet #' + str(counter) + ' of ' + str(size) + ' *****') 
        print('OLD TRIPLET:')
        print(tup[0].toStr())
        newTriplet = addDetections(tup[0], tup[1], expDict)
        print('NEW TRIPLET:')
        print(newTriplet.toStr())
        newTripList.append(newTriplet)
    name = 'crossCampaignTriplets+' + saveName
    writeTriplets(newTripList, name + '.txt', True)
    pickleTriplets(newTripList, name + '.pickle')
'''
if __name__=='__main__':
    main()
