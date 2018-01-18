import pickle
import argparse
import time
import sys
sys.path.insert(0, './linkerCode')
from LinkerLib import Triplet
from LinkerLib import pickleTriplets
from LinkerLib import writeTriplets
from LinkerLib import printPercentage
#splits a list of triplets into 'chunks' number of chunks
#plus one extra if #chunks doesn't divide the size of the triplets
'''
input: list of triplets, number of chunks to split the triplets
output: a list of list of triplets
'''
def splitList(tripList, chunks):
    print('spliting list of size ' + str(len(tripList)) + 
        ' into ' + str(chunks) + ' sized chunks')
    size = len(tripList)
    #sorry i made chunks be the chunksize instead of 
    # number of chunks i promise i'll fix the variables
    # and documentation soon
    chunkSize = chunks
    chunks = size/chunkSize
    print('chunk size: ' + str(chunkSize))
    splitTrips = []
    time0 = time.time()
    for x in range(chunks):
        chunkList = []
        printPercentage(x, chunks, time.time()-time0)
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
    args.add_argument('triplets', nargs=1, help='path to pickle file; ' + 
             'file has format triplet+SNOBS_SEASON###_ML0#.pickle')
    args.add_argument('numChunks', nargs='?', default=100, 
        help='number of chunks')
    args = args.parse_args()
    print('open pickle file')
    pickleFile = open(args.triplets[0], 'rb')
    print('loading pickle file')
    triplets = pickle.load(pickleFile)
    tripChunks = splitList(triplets, int(args.numChunks))
    
    #save triplets
    saveName = args.triplets[0].split('+')[-1].split('.')[0]
    saveName2 = args.triplets[0].split('+')[-2]
    if '/' in saveName2:
        saveName2 = saveName2.split('/')[-1]
    print(saveName2)
    for x in range(len(tripChunks)):
        writeTriplets(tripChunks[x], 'chunk' + str(x) + '+' + saveName2 + '+' + saveName + '.txt', False)
        pickleTriplets(tripChunks[x], 'chunk' + str(x) + '+' + saveName2 + '+' + saveName + '.pickle')  

if __name__ == '__main__':
    main() 
    
