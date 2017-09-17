import sys
sys.path.insert(0, './linkerCode/')

import argparse
import pandas as pd
import time
from LinkerLib import printPercentage
# removes elements with ccd
def rmCCDs(df, thresh):
    print('number of detections: ' + str(len(df)))
    #a dictionary of dictionaries
    ccdict = {}
    expnums = df['EXPNUM'].tolist()
    ccdnums = df['CCDNUM'].tolist()
    time0 = time.time()
    #count the number of things with certain ccd
    for x in range(len(df)):
        printPercentage(x, len(df), time.time()-time0)
        exp = expnums[x]
        ccd = ccdnums[x]
        if exp in ccdict:
            if ccd in ccdict[exp]:
                ccdict[exp][ccd] += 1
            else:
                ccdict[exp][ccd] = 1 
        else:
            ccdict[exp] = {}
    #list of (expnum, ccdnum)
    overThresh = []
    #check if any are over thresh
    for key, value in ccdict.iteritems():
        for key2, value2 in value.iteritems(): 
            if value2 > thresh:
                print('expnum:' + str(key) + ' ccdnum:' + str(key2) + ' hits:' + str(value2))
                overThresh.append((key, key2))
    #to do remove detections with corresponding overthresh
    rmList = []
    for x in range(len(df)):
        if ((expnums[x], ccdnums[x]) in overThresh):
            rmList.append(x)
    print('rows removed:' + str(len(rmList)))
    df = df.drop(df.index[rmList])
    print('size after remove: ' + str(len(df)))
    return df

def main():
    args = argparse.ArgumentParser()
    #pass in varRM+ file
    args.add_argument('detections', nargs=1, help='path to csv detection file; ' + 
                        'filename has format varRM+SNOBS_SEASON###_ML0#.csv')
    args.add_argument('ccdThresh', nargs='?', default=100, help='threshold for disregarding ccd')
    args = args.parse_args()
    savename = args.detections[0].split('+')[-1].split('/')[-1].split('.')[0]
    df = pd.read_csv(args.detections[0])
    df = rmCCDs(df, int(args.ccdThresh))
    df.to_csv('ccdRM+' + savename + '.csv')

if __name__=='__main__':
    main()
