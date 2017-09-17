'''
under construction
attempts to take lists of triplets and merge them into longer orbits

'''
import time
from LinkerLib import Triplet
from LinkerLib import Detection
from LinkerLib import printPercentage
from LinkerLib import writeTriplets
from LinkerLib import pickleTriplets
from GraphNodes import connectGraph
from GraphNodes import graphLinks
import argparse
import pickle
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
import itertools
   
#removes detections one by one to test if orbits fit
def removeDets(triplets):
    goodTriplets = []
    badTriplets = []
    time0 = time.time()
    counter = 0
    for trip in triplets:
        printPercentage(counter, len(triplets), time.time()-time0)
        elements, errs= trip.calcOrbit()
        if(elements['a'] > 2 and elements['e'] < 1):
            goodTriplets.append(trip)
        elif(len(trip.dets) > 3):
            found = False
            for x in range(len(trip.dets)):
                tripTest = Triplet(trip.dets[:x] + trip.dets[(x+1):])
                elements, errs = tripTest.calcOrbit()
                if(elements['a'] > 2 and elements['e'] < 1):
                    goodTriplets.append(tripTest)
                    found = True
                    continue

            if(not found):
                badTriplets.append(trip)
            
        else:
            badTriplets.append(trip)
        counter+=1
    return goodTriplets, badTriplets            

#given a list of triplets, return a list of triplets where no two detections have the same exposure
def filterSameExp(triplets):
    print('\nfiltering same exposures')
    result = []
    time0 = time.time()
    count = 0
    for trip in triplets:
        count+=1
        printPercentage(count, len(triplets), time.time()-time0)
        dets = trip.dets
        #place each detection in its respective night
        mjds = set([x.mjd for x in dets]) 
        #dictionary where mjd is key, and value is list of detections
        mjdCount = dict([(x,[]) for x in mjds])
        for det in dets:
            mjdCount[det.mjd].append(det)
        #convert dictionary to list
        countList = mjdCount.values()
        #get every combination of those triplets
        result += (list(itertools.product(*countList)))
    return result

#returns list of triplets with good paramenters
def potentialTriplets(triplets):
    size = len(triplets)
    goodList = []
    badList = []
    time0 = time.time()
    counter = 0
    for trip in triplets:
        counter += 1
        printPercentage(counter, size, time.time()-time0)
        elements, errs = trip.calcOrbit()
        if(elements['a'] > 2 and elements['e']<1):
            goodList.append(trip)
        else:
            badList.append(trip)
    return goodList, badList

#given a graph, finds the connected components and attempts to identify objects from them
def findObjects(graph):
    print('getting connected components')
    components = nx.connected_components(graph)#.to_undirected())
    print('wrapping to triplet')
    triplets = [Triplet(list(comp)) for comp in components]
    print('done wrapping')
    '''
    #remove same exposure detections
    diffExps = filterSameExp(triplets)
    print('wrapping more triplets')
    diffExptrips = [Triplet(diffExp) for diffExp in diffExps]
    print(len(diffExptrips))
    '''
    #test remove one or two detections
    print('check removing detections')
    goodTriplets, badTriplets = potentialTriplets(triplets)#removeDets(triplets)#diffExptrips)
    return goodTriplets, badTriplets        
   
def main():
    args = argparse.ArgumentParser()
    args.add_argument('triplets', nargs=1, help='path to pickle file')
    args.add_argument('graph', nargs='?', default='n', help='whether to graph, [y] for yes')
    args = args.parse_args()
    print('loading pickle file')
    triplets = pickle.load(open(args.triplets[0],'rb'))
    saveName = args.triplets[0].split('+')[-1].split('.')[0]
    # merges all triplets onto a graph
    print('placing on graph')
    G, dict = connectGraph(triplets)
    
    print('\nfinding objects in graphs')
    #extracts orbits out of merged triplets
    good, bad = findObjects(G)

    #saves every merged triplets 
    saveGoodName = 'goodMergedTriplets+' + saveName
    saveBadName = 'badMergedTriplets+' + saveName 
    writeTriplets(good, saveGoodName + '.txt', True)
    writeTriplets(bad, saveBadName + '.txt', True)
    pickleTriplets(good, saveGoodName + '.pickle')
    pickleTriplets(bad, saveBadName + '.pickle')
    
    if(args.graph == 'y'): 
        #graphs original merged graph
        graphLinks(G, dict, saveName)#, [-20, 0,-52,-40])
        G, dict= connectGraph(good)
        graphLinks(G, dict, saveName)
    #, [3.4,5.2,-1,0])
    '''
    merged = findObjects(G)
    
    savename = 'mergedTrips+' + saveName + '.txt'
    writeTriplets(merged, savename)
    
    print('writing pickle file')
    savename = 'mergedTrips+' + saveName + '.pickle'
    with open(savename, 'wb') as f:
        pickle.dump(merged,f)
    '''
if __name__ == '__main__':
    main()
