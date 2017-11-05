import sys
sys.path.insert(0, 'linkerCode/')
import pickle
import argparse
import time
from combineNLet import getBounds, partition, combine, removeDuplicates, removeSameExposure, checkGoodOrbit
from LinkerLib import Triplet, writeTriplets, pickleTriplets

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('nlets', nargs=1, help='path to nlet pickle file')
    parser.add_argument('chisq', nargs='?', default = 5,
                            help='chisq threshold for good orbit')
    args = parser.parse_args()
    
    nlets = pickle.load(open(args.nlets[0], 'rb'))
    savename = args.nlets[0].split('+')[-1].split('.')[0]
    
    print('Recursively merging...')

    while len(nlets) != 0:
        t0 = time.time()
        numdet = len(nlets[0].dets)
        print('Processing ' + str(len(nlets)) + ' ' + str(numdet) +'-lets')
        minRA, maxRA, minDec, maxDec = getBounds(nlets)
        arr = partition(nlets, minRA, maxRA, minDec, maxDec, 700, 700)
        nPlusOne, unchanged = combine(arr, numdet)
        
        # remove duplicates
        nPlusOne = removeDuplicates(nPlusOne)
        unchanged = removeDuplicates(unchanged)

        # remove same exposure
        nPlusOne = removeSameExposure(nPlusOne)
        unchange = removeSameExposure(unchanged)

        # check for good orbits
        nPlusOne = checkGoodOrbit(nPlusOne, float(args.chisq))
        unchanged = checkGoodOrbit(unchanged, float(args.chisq))

        num_u = len(unchanged)

        if num_u != 0:
            print('Saving ' + str(num_u) + ' unchanged ' + str(numdet) + '-lets')
            saveas_u = 'good' + str(numdet) + 'letU' + str(args.chisq) + '+' + savename
            writeTriplets(unchanged, saveas_u + '.txt', False)
            pickleTriplets(unchanged, saveas_u + '.pickle')
            
        nlets = nPlusOne
        print('Time taken this cycle: ' + str(time.time() - t0))
    print('Done')

if __name__ == '__main__':
    main()

