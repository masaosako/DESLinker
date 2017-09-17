import sys
sys.path.insert(0, '../csvDetectionfiles')
sys.path.insert(0, '../')

import numpy as np
from LinkerLib import Detection

#takes two detections and determines if their mags are within range
#subject to change with more incoming data

def withinMag(det1, det2):
    #very lenient cut
    magThresh = 1.0
    if(det1.band == det2.band):
        return abs(det1.mag-det2.mag) < magThresh
    else:
        bandDict = {'z': 0, 'i': 0.375, 'r': 0.875, 'g': 1.75}
        mag1 = det1.mag - bandDict[det1.band]
        mag2 = det2.mag - bandDict[det2.band] 
        #print('mag1: ' + str(mag1) + ', mag2: ' + str(mag2))
        return abs(mag1-mag2) < magThresh

def withinSpeed(det1, det2, det3):
    #calculating speed
    deltaRa1 = det2.ra - det1.ra
    deltaDec1 = det2.dec - det1.dec
    deltaRa2 = det3.ra - det2.ra
    deltaDec2 = det3.dec - det2.dec
    dist1 = np.sqrt(deltaRa1*deltaRa1 + deltaDec1*deltaDec1)
    dist2 = np.sqrt(deltaRa2*deltaRa2 + deltaDec2*deltaDec2)
    deltaMJD1 = det2.mjd - det1.mjd
    deltaMJD2 = det3.mjd - det2.mjd
    speed1 = dist1/deltaMJD1
    speed2 = dist2/deltaMJD2
    #print("speed1/speed2: " + str(speed1/speed2))
    #print("speed2/speed1: " + str(speed2/speed1))
    #ratio of speeds within speedThresh
    # Originally 0.15 + 0.02*abs(deltaMJD1 + deltaMJD2)
    speedThresh = 0.15 + 0.02*abs(deltaMJD1+deltaMJD2)
    #print('speedThresh: ' + str(speedThresh))
    #print('Delta MJD: ' + str(deltaMJD1 + deltaMJD2))
    #print('\n')
    return speed1/speed2 < speedThresh+1 and speed2/speed1 < speedThresh+1

def withinInc(det1, det2, det3):
    return True
