import sys

import time
from LinkerLib import Detection
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
 
#form a list of triplets given a list of detections with links
def formTriplets(args):
    detPairs, queue = args
    #queue = 0 if not multiprocessing
    tripList = []
    time0 = time.time()
    counter = 0
    for det in detPairs:
        if(queue == 0):
            counter += 1
            printPercentage(counter, len(detPairs), time.time()-time0)
        for link in det.linkedList:
            for trip in link.linkedList:
                triplet = Triplet([det, link, trip])
                tripList.append(triplet)
    if(queue != 0):
        queue.put(tripList)
    return tripList

#turns pairs into triplets by multiprocessing
'''
input: list of detections with their links
        number of cpus to be used
output: list of triplets formed by the links
'''
def multiProcessPairs(detPairs, Ncpu):
    print('forming triplets')
    pool = Pool(Ncpu)
    queue = Manager().Queue()
    splitList = [([det], queue) for det in detPairs]
    result = pool.map_async(formTriplets,splitList)
    time0 = time.time()
    
    #monitor the loop 
    while not result.ready():
        print('triplet forming progress: ' + str(queue.qsize()) + 
            ' of ' + str(len(detPairs)) + ' triplets form ' + 
                str(float(queue.qsize())/len(detPairs)*100) + '% at ' + 
                str(time.time()-time0) + ' sec')
        time.sleep(1)
    print('ready: ' + str(result.ready()) + ' size: ' + str(queue.qsize()))
    
    triplets = []
    real_result = result.get()
    #compile result into a list
    for r in real_result:
        triplets.append(r) 
    
    triplets = [trip for trips in triplets for trip in trips]
    return triplets

#splits a list of triplets into 'chunks' number of chunks
#plus one extra if #chunks doesn't divide the size of the triplets
'''
input: list of triplets, number of chunks to split the triplets
output: a list of list of triplets
'''
def splitList(tripList, chunks):
    print('\nspliting list of size ' + str(len(tripList)) +
        ' into ' + str(chunks) + ' chunks')
    size = len(tripList)
    #rounded down
    chunkSize = int(size/int(chunks))
    splitTrips = []
    time0 = time.time()
    for x in range(int(chunks)):
        chunkList = []
        printPercentage(x, int(chunks), time.time()-time0)
        for y in range(chunkSize):
            chunkList.append(tripList.pop())
        splitTrips.append(chunkList)
    print('\nappending ' + str(len(tripList)) + ' remaining trips')
    #the remainder of size/chunks
    if(tripList):
        splitTrips.append(tripList)
    return splitTrips


def main():
    args = argparse.ArgumentParser()
    args.add_argument('linkedPairs', nargs=1, 
                        help='path to pickle file; ' +
                        'file has format detectionLinks+SNOBS_SEASON###_ML0#.pickle')
    args.add_argument('numChunks', nargs='?', default=1, help='number of chunks')
    args.add_argument('ncpus', nargs='?', default=1, 
                        help='number of cpus')
    args = args.parse_args()
    
    #load the pickle file with list of detections and links
    print('load pickle file')
    detPairs = pickle.load(open(args.linkedPairs[0], 'rb'))
    
    #makes a list of potential triplets
    Ncpu = int(args.ncpus)
    triplets = formTriplets((detPairs, 0))
    numChunks = args.numChunks
    tripChunks = splitList(triplets, numChunks)
    #triplets = multiProcessPairs(detPairs, Ncpu)    
    saveName = args.linkedPairs[0].split('+')[-1].split('.')[0]
    for x in range(len(tripChunks)):
        writeTriplets(tripChunks[x], 'chunk{0:04d}'.format(x) + '+' + saveName + '.txt', False)
        pickleTriplets(tripChunks[x], 'chunk{0:04d}'.format(x) + '+' + saveName + '.pickle')


if __name__ == '__main__':
    main()


