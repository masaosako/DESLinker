'''
under construction
a second version of mergeTrips.py that attempts to merge triplets into longer orbits
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

def contains(small, large):
    for ele in small:
        if ele not in large:
            return False
    return True

def worthadding(newblob, listofblob):
    count = 0
    for oldblob in listofblob:
        oldblob = oldblob.dets
        if contains(newblob.dets, oldblob):
            return listofblob # oldblob already has it
        if contains(oldblob, newblob.dets):
            # should put the new blob in and old blob out
            listofblob.pop(count)
            listofblob.append(newblob)
            return listofblob
        count = count + 1
    listofblob.append(newblob)
    return listofblob

def growtrees(graph, d):
    good = []
    for CC in d:
        CCgood = []
        for trip in d[CC]:
            blobsnow = [trip]
            while blobsnow:
                for blob in blobsnow:
                    successors = graph.successors(blob.dets[-1])
                    if not successors:
                        blobsnow.remove(blob)
                        worthadding(blob, CCgood)
                        continue
                    for next in successors:
                        trydets = blob.dets[:]
                        trydets.append(next)
                        tryblob = Triplet(trydets)
                        elements, errs = tryblob.calcOrbit()
                        if(elements['a'] > 2 and elements['e'] < 1):
                            if blob in blobsnow:
                                blobsnow.remove(blob)
                            blobsnow.append(tryblob)
                    worthadding(blob,CCgood)
                    if blob in blobsnow:
                        blobsnow.remove(blob)
        good.extend(CCgood)
    return good

#given a graph, finds the connected components and attempts to identify objects from them
def organizetriplets(graph, triplets):
    components = list(nx.connected_component_subgraphs(graph.to_undirected(),copy=False))
    # This dictionary associates triplets with the connect components they belong to
    d = dict([(com,[]) for com in components])
    for trip in triplets:
        for com in components:
            if com.has_node(trip.dets[0]):
                #print(com, ' has ', trip, '\n')
                d[com].append(trip)
                break
    return d        
   
def main():
    args = argparse.ArgumentParser()
    args.add_argument('triplets', nargs=1, help='path to pickle file')
    args = args.parse_args()
    triplets = pickle.load(open(args.triplets[0],'rb'))
    saveName = args.triplets[0].split('+')[-1].split('.')[0]
    # merges all triplets onto a graph
    G, dict = connectGraph(triplets)
    
    # know the tripelts inside each connected components
    d  = organizetriplets(G, triplets)

    good = growtrees(G,d)

    #saves every merged triplets 
    saveGoodName = 'goodMergedTriplets3+' + saveName
    #saveBadName = 'badTriplets+' + saveName 
    writeTriplets(good, saveGoodName + '.txt')
    #writeTriplets(bad, saveBadName + '.txt')
    pickleTriplets(good, saveGoodName + '.pickle')
    #pickleTriplets(bad, saveBadName + '.pickle')
    
    # graphs original merged graph
    #graphLinks(G, dict, saveName)
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
