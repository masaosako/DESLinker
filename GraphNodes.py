'''
The library related to any part of the pipeline that requires graphing the data
and displaying that graph

'''




import sys
#sys.path.insert(0,'../')
sys.path.insert(0,'../csvDetectionFiles')

import pickle
import time
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
from LinkerLib import printPercentage
from LinkerLib import Detection
from LinkerLib import Triplet

font = {'family' : 'normal',
            'weight' : 'bold',
            'size' : 22}
matplotlib.rc('font', **font)

#graphs linkpairs
'''
input: detections to be linked, whether or not to add edges to graph
output: graph G and positions dict
'''
def linkGraph(detections, addEdges=True):
    G = nx.Graph()
    startT = time.time()
    count = 0
    dict = {}
    for det in detections:
        printPercentage(count, len(detections), time.time()-startT)
        count += 1
        G.add_node(det)
        dict[det] = (det.ra, det.dec)
        for link in det.linkedList:
            G.add_node(link)
            dict[link] = (link.ra, link.dec)
            if(addEdges):
                G.add_edge(det, link)
    return G, dict


# forms a graph of triplets
'''
input: list of Triplets
output: a graph G containing nodes and edges
        a dictionary dict containing positions for graphing
'''
def connectGraph(triplets):
    G = nx.Graph()
    startT = time.time()
    count = 0
    dict = {}
    for trip in triplets:
        printPercentage(count, len(triplets), time.time()-startT)
        count += 1
        detList = sorted(trip.dets, key=lambda detection: detection.mjd)
        G.add_node(detList[0])
        dict[(detList[0])] = (detList[0].ra, detList[0].dec)
        for i in range(1,len(detList)):
            G.add_node(detList[i])
            dict[(detList[i])] = (detList[i].ra, detList[i].dec)
            G.add_edge(detList[i-1], detList[i])
    return G, dict

#takes a graph G, positions in dict, and saves to saveName
#optional axes
#lays out the graph
'''
input: graph G with nodes of detections and edges of links
        dict with position of nodes
        saveName: title of graph
        axes: optional, boundary of the graph
output: saves and shows the graph
'''
def graphLinks(G, dict, saveName, axes=[0,0,0,0]):
    print('\nlaying out graph')
    fixed_nodes = dict.keys()
    print('found keys')
    time2 = time.time()
    pos = nx.spring_layout(G, pos=dict, fixed=fixed_nodes)
    time3 = time.time()
    print('done layout after ' + str(time3-time2) + 'sec')
    print('saving positions as pickle file')
    with open('graphPos.pickle', 'wb') as f:
        pickle.dump((pos,G),f)
    print('saved')
    #show and save the graph
    graphPos(pos, G, saveName, axes)

#graphs G with positions pos
#does the actual graphing
'''
input: graph G with nodes of detections and edges of links
        pos with position of nodes
        axes: optional, boundy of graph
        saveName, boundry of graph
output: saves and shows the graph
'''
def graphPos(pos, G, saveName, axes=[0,0,0,0]):
    print('setting figsize')
    plt.figure(figsize=(12,12))
    print('checking axes')
    if(axes!=[0,0,0,0]):
        plt.axis(axes)
    print('drawing network')
    time4 = time.time()
    nx.draw_networkx(G, pos, node_size=10, font_size=8)
    print('done drawing after' + str(time.time()-time4) + 'sec')
    plt.title(saveName)
    plt.xlabel('RA (degrees)')
    plt.ylabel('Dec (degrees)')
    print('number of nodes: ' + str(len(G)))
    graphName = 'merged+' + saveName + '.svg'
    plt.savefig(graphName, format='svg', dpi=1200)
    plt.show()

#graphs things if graphing terminates before it's finished
def main():
    args = argparse.ArgumentParser()
    args.add_argument('graph', nargs=1, help='path to pickle file')
    args = args.parse_args()
    print('opening pickle file: ' + args.graph[0])
    (pos, G) = pickle.load(open(args.graph[0],'rb'))
    saveName = args.graph[0].split('+')[-1].split('.')[0]
    print('pickle file loaded')
    graphPos(pos, G, saveName)#, [3.4,5.2,-1,0])
    
if __name__=='__main__':
    main() 
