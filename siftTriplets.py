import sys
sys.path.insert(0, './linkerCode')
import time
from LinkerLib import printPercentage
from LinkerLib import Triplet
from LinkerLib import writeTriplets
from LinkerLib import pickleTriplets

import numpy as np
import pickle
import argparse
import time
import ephem

from multiprocessing import Pool, cpu_count, Manager, Queue

#find the potential orbits given a list of tripletsi
'''
input: a list of triplets
output: goodList, a list of triplets that satisfy orbits
        badList, a list of triplets that don't satisfy orbits
        missList, a list of triplets that don't satisfy, but are P9 like or
            are the same fake for fakes
'''
def potentialTriplets(args):#, queue):
    triplets, queue= args
    #queue is 0 if no multiprocessing
    size = len(triplets)
    counter = 0
    time0 = time.time()
    goodList = []
    badList = []
    p9List = []
    missList = []
    # for each detection, loop through each possible combination of triplets
    for trip in triplets:
       # if(queue == 0):
            #printPercentage(counter, size, time.time()-time0)

        elements, errs = trip.calcOrbit()
        #chisq = trip.getChiSq()
        #orbit is a swigpy object that can't be saved for some reason
        trip.orbit = 0
        if(elements['a'] > 2 and elements['e'] < 1):
        #if chisq < threshold:
            goodList.append(trip)
        else:
            # already missing quite a lot before this stage
            if(trip.p9Like()):
                p9List.append(trip)
            elif(trip.sameFake()):
                missList.append(trip)
            else:
                badList.append(trip)
        counter += 1
    badid = []
    for bad in badList:
        for det in bad.dets:
            badid.append(det.fakeid)

    missid = [x.dets[0].fakeid for x in missList]
    goodid = [x.dets[0].fakeid for x in goodList]
    p9id = [x.dets[0].fakeid for x in p9List]
    #print('\n Number missing: ' + str(len(set(missid))))
    #print('Number good: ' + str(len(set(goodid))))
    #print('Number p9 like: ' + str(len(set(p9id))))
    #print('Number bad: ' + str(len(set(badid))))
    #print('Total: ' + str(len(set(badid+missid+goodid))))i
    if(queue != 0):
        queue.put(allTrips)
    return goodList, missList, badList, p9List

#check triplets with multiproccessing
def multiProcessTriplets(triplets, Ncpu):
    print('sifting triplets')
    pool = Pool(Ncpu)
    manager = Manager()
    queue = manager.Queue()
    splitList = [([triplet],queue) for triplet in triplets]
    result = pool.map_async(potentialTriplets, splitList)
    #monitor loop
    size_old = 0
    time0 = time.time()
    while not result.ready():
        print('triplet filtering progress: ' + str(queue.qsize()) +
            ' of ' + str(len(triplets)) + ' triplets checked ' +
                str(float(queue.qsize())/len(triplets)*100) + '% at ' +
                str(time.time()-time0) + ' sec')
        time.sleep(1)
    
    #consolidate reults into a few lists
    goodTrips = []
    missTrips = []
    badTrips = []
    p9Trips = []
    allTrips = []
    for r in result.get():
        goodTrip, missTrip, badTrip, p9Trip, allTrip = r
        goodTrips.append(goodTrip)
        missTrips.append(missTrip)
        badTrips.append(badTrip)
        p9Trips.append(p9Trip)
        allTrips.append(allTrip)
    pool.close()
    goodTrips = [trip for triplets in goodTrips for trip in triplets]
    missTrips = [trip for triplets in missTrips for trip in triplets]
    badTrips = [trip for triplets in badTrips for trip in triplets]
    p9Trips = [trip for triplets in p9Trips for trip in triplets]
    allTrips = [trip for triplets in allTrips for trip in triplets]
    return goodTrips, missTrips, badTrips, p9Trips, allTrips

def main():
    args = argparse.ArgumentParser()
    args.add_argument('triplets', nargs=1,
                        help='path to pickle file; file has format ' + 
                        'chunk###+SNOBS_SEASON###_ML0#.pickle')
    #args.add_argument('chisq', nargs='?', default=5,
    #                    help='chisq threshold to be considered a good triplet')
    args.add_argument('ncpus', nargs='?', default=1, 
                        help='number of cpus')
    args = args.parse_args()
    print('open pickle file ' + args.triplets[0]) 
    pickleFile = open(args.triplets[0], 'rb')
    print('load pickle file')
    triplets = pickle.load(pickleFile)
    Ncpu = int(args.ncpus)

    #checks orbits for good fit
    goodTrips, missTrips, badTrips, p9Trips = potentialTriplets((triplets, 0))
    #multiProcessTriplets(triplets, Ncpu)


    #saving triplets
    saveName = args.triplets[0].split('+')[-1].split('.')[0]
    chunkName = args.triplets[0].split('/')[-1].split('+')[0]
    saveName = chunkName + '+' + saveName
    writeTriplets(goodTrips, 'goodtriplets+' + saveName + '.txt', False)
    writeTriplets(missTrips, 'misstriplets+' + saveName + '.txt', False)
    writeTriplets(badTrips, 'badtriplets+' + saveName + '.txt', False)
    writeTriplets(p9Trips, 'p9triplets+' + saveName + '.txt', False)

    pickleTriplets(goodTrips, 'goodtriplets+' + saveName + '.pickle')
    pickleTriplets(missTrips, 'misstriplets+' + saveName + '.pickle')
    pickleTriplets(badTrips, 'badtriplets+' + saveName + '.pickle')
    pickleTriplets(p9Trips, 'p9triplets+' + saveName + '.pickle')

if __name__=='__main__':
    main()
