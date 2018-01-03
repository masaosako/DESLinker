import sys
sys.path.insert(0,'linkerCode/')
import pickle
import argparse
import LinkerLib
from LinkerLib import Detection
from LinkerLib import printPercentage
import time
import pandas as pd

def getYear(mjd):
    if 56505 <= mjd <= 56717:
        return 0
    elif 56870 <= mjd <= 57082:
        return 1
    elif 57235 <= mjd <= 57448:
        return 2
    elif 57601 <= mjd <= 57813:
        return 3
    print("ERROR: nlet mjd not within campaign bounds:" + str(mjd))
    return -1

# assumes nlets from single season
def detPrediction(nlets, filename):
    df = pd.read_csv(filename, header=None, delimiter=',',
                    names=['ind', 'expnum', 'nite', 'mjd', 'ra', 'dec',
                            'band', 'exptime', 'on', 'num1', 'num2', 'num3'])
    explist = df['expnum']
    #ralist = df['ra']
    #declist = df['dec']
    mjdlist = df['mjd']
    
    year = getYear(nlets[0].dets[0].mjd)

    division = {0:[], 1:[], 2:[], 3:[]}

    for i in range(len(explist)):
        index = getYear(mjdlist[i])
        if index == -1:
            continue
        division[index].append((explist[i], mjdlist[i]))

    if year - 1 in division:
        prev = division[year - 1]
    else:
        prev = []
    if year + 1 in division:
        following = division[year + 1]
    else:
        following = []

    searchSpan = prev + following + division[year]

    predDets = []
    counter = 0
    nletsSize = len(nlets)
    for nlet in nlets:
        counter += 1
        predictions = []
        count = 0
        time0 = time.time()
        size = len(searchSpan)
        print('\n triplet number ' + str(counter) + ' of ' + str(nletsSize) + ' predictions for:')
        print(nlet.toStr())
        for i in range(size):
            count += 1
            printPercentage(count, size, time.time()-time0)
            coord, erra, errb, pa = nlet.predictPos(searchSpan[i][1])
            det = Detection(coord[0], coord[1], searchSpan[i][1],
                              2, 0, searchSpan[i][0], 0, 0, 60, 0)
            det.erra = erra
            det.errb = errb
            det.pa = pa
            predictions.append(det)
        predDets.append(predictions)
    return zip(nlets, predDets)

    
        

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('nlets', nargs=1, help='path to nlets')
    parser.add_argument('listfile', nargs=1, help='path to detection link file')
    args = parser.parse_args()

    print('Loading nlets')
    nlets = pickle.load(open(args.nlets[0], 'rb'))
    print('Nlets loaded')

    #df = pd.read_csv(args.listfile[0], header=None, delim_whitespace=True,
    #                    names=['ind', 'expnum', 'nite', 'mjd', 'ra', 'dec',
    #                            'band', 'exptime', 'on', 'num1', 'num2', 'num3'])
    print('Forming predicted Detection objects')
    predDets = detPrediction(nlets[:10], args.listfile[0])
    savename = args.nlets[0].split('+')[-1].split('.')[0]
    savename2 = args.nlets[0].split('+')[-2].split('/')[-1]
    savename='listOfPredictions+' + savename2 + '+' + savename
    print('\nsaving predictions')
    with open(savename + '.pickle', 'wb') as f:
        for x in predDets:
            x[0].orbit=0
        pickle.dump(predDets, f)
    '''
    for predDet in predDets:
        print(predDet[0].toStr())
        print('***************************')
        for det in predDet[1]:
            print(det.toStr())
    '''

    '''
    detlink = pickle.load(open(args.detlink[0], 'rb'))
    expdict = expDictionary(detlink)

    for k,v in expdict.items():
        print("Expnum: " + str(k))
        print("---------------------")
        for det in v:
            print(det.toStr())
        print('')
    '''


if __name__ == '__main__':
    main()
