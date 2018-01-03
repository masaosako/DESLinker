from LinkerLib import expDictionary
from LinkerLib import Triplet
from LinkerLib import printPercentage
from LinkerLib import pickleTriplets, writeTriplets
from generateDetections import wrapDets
from fitMoreDets import addToNlet
import numpy as np
import time
import pickle
import argparse


def main():
    args = argparse.ArgumentParser()
    args.add_argument('detList', nargs=1, help='path to csv list of detections')
    args.add_argument('potentialNlets', nargs=1, help='Nlets we want to add on to')
    args = args.parse_args()
    print('\nloading nlets ' + args.potentialNlets[0])
    nlets = pickle.load(open(args.potentialNlets[0], 'rb'))
    saveName = args.potentialNlets[0].split('+')[-1].split('.')[0]
    saveName2 = args.potentialNlets[0].split('+')[-2]
    if '/' in saveName2:
        saveName2 = saveName2.split('/')[-1]
    print(saveName2)

    print('\nwrapping detections: ' + args.detList[0])
    detList = wrapDets(args.detList[0])
    print('\nmaking dictionary')
    expDict = expDictionary(detList)
    size = len(nlets)
    counter = 0
    time0 = time.time()
    newNletList = []
    for nlet in nlets:
        counter += 1
        print('\n*****adding to Nlet #' + str(counter) + ' of ' + str(size) + ' *****')
        print('OLD TRIPLET:')
        print(nlet.toStr())
        newLet = addToNlet(nlet, expDict)
        print('\nNEW TRIPLET:')
        print(newLet.toStr())
        newNletList.append(newLet)
    name = 'crossCampaignTriplets+' + saveName2 + '+' + saveName
    writeTriplets(newNletList, name + '.txt', True)
    pickleTriplets(newNletList,name + '.pickle')
    
if __name__ == '__main__':
    main()
