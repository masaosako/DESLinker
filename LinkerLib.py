import sys
sys.path.append('..')
import numpy as np
import GammaTPlotwStatTNOExFaster as gts
import ephem
from astropy.time import Time
from datetime import datetime
import sys
from Orbit import Orbit
import time
import pickle
from orbfitScript import degToHour, degToDMS

class Region:
    # fields include the range for RA / Dec and list of detections
    def __init__(self, raRng, decRng):
        self.raLo = raRng[0]
        self.raHi = raRng[1]
        self.decLo = decRng[0]
        self.decHi = decRng[1]
        self.detections = []
        
    def add(self, tno):
        self.detections.append(tno)
    
    def toStr(self):
        minbounds = '(minRa: ' + str(self.raLo) + ', minDec: ' + str(self.decLo) + ')'
        maxbounds = '(maxRa: ' + str(self.raHi) + ', maxDec: ' + str(self.decHi) + ')'
        return minbounds + '  ' + maxbounds 

#all the parameters of a detection
#ra is between -180 and 180
class Detection:
    def __init__(self, ra, dec, mjd, flux, objid, expnum, ccd, band, lookAhead, fakeid):
        self.ra = ra
        if self.ra > 180:
            self.ra = ra -360
        self.dec = dec
        self.mjd = mjd
        self.mag = -2.5*np.log10(flux) + 31.4
        self.objid = objid
        self.expnum = expnum
        self.ccd = ccd
        self.band = band
        self.fakeid = fakeid
        self.lookAhead = lookAhead
        #describes a cone region that would contain possible pair links
        #angle in radians where 0 is along the x-axis, where change in dec is 0
        #angle1 < angle2
        # is an array of angles and radii for each nite with index corresponding
        # to the nite ahead
        self.angle1, self.angle2, self.radius = gts.calcCone(lookAhead, 
                        self.ra, self.dec, self.mjd)
        #list of linked detections (intially empty)
        self.linkedList = []
    #so nodes can have class of Detection
    def __hash__(self):
        return hash(self.objid)
    def __eq__(self, other):
        return self.objid == other.objid
    def __str__(self):
        return ''#str(self.fakeid%1000)
    def toStr(self):
        coord = '(RA: ' + str(self.ra) + ', DEC: ' + str(self.dec) + ') '
        id = ('fakeid: ' + str(self.fakeid) + ', objid: ' + str(self.objid) + 
                ', expnum: ' + str(self.expnum) + ', ccd: ' + str(self.ccd))
        return (id + ', mjd: ' + str(self.mjd) + 
                ', mag: ' + str(self.mag) + ', band: ' + str(self.band) 
                + ' ' + coord)
    def toDat(self):
        # formerly used str(ephem.hours(self.ra*np.pi/180))
        # str(ephem.degrees(self.dec*np.pi/180))
        s = (str(self.mjd+2400000.5) + ' ' + degToHour(self.ra) +
                ' ' + degToDMS(self.dec) + ' 0.1 807')
        return s
    #checks if another detection is inside the cone
    def withinCone(self, det):
        nitesAhead = int(det.mjd - self.mjd)
        #only looking ahead a certain number of nites
        if(nitesAhead + 1 >= self.lookAhead or nitesAhead < 0):
            return False
        ra = det.ra
        # make sure ra is between -180 and 180
        if ra > 180:
            ra = ra - 360
        deltaRa = ra-self.ra
        deltaDec = det.dec-self.dec
        dist = np.sqrt(deltaRa**2 + deltaDec**2)
        angle = np.arctan2(deltaDec, deltaRa)
        #vector angle is between 0 and 360
        if(angle < 0):
            angle += 2*np.pi
        #checks if the cone overlaps with the 0 angle
        rev = self.angle1[nitesAhead] > self.angle2[nitesAhead]
        #a one night interval
        angle1 = max(self.angle1[nitesAhead], self.angle1[nitesAhead + 1])
        angle2 = min(self.angle2[nitesAhead], self.angle2[nitesAhead + 1])
        radius = max(self.radius[nitesAhead], self.radius[nitesAhead + 1])
        if(rev):
            if(angle < angle1 and angle > angle2 and dist < radius):
#                print('\nself: ' + self.toStr() + '\n  them: '  + det.toStr())
#                print('dist: ' + str(dist))
                return True
        
            else:
                return False
        #the cone overlaps with the 0 angle
        else:
            if((angle < angle1 or angle > angle2) and dist < radius):
                return True
            else:
                return False

    #returns the min/max ra/dec of the cone
    def bounds(self): 
        vertex1 = [self.ra, self.dec]
        vertex2 = [self.ra + np.amax(self.radius)*np.cos(np.amax(self.angle1)),
                    self.dec + np.amax(self.radius)*np.sin(np.amax(self.angle1))]
        vertex3 = [self.ra + np.amax(self.radius)*np.cos(np.amin(self.angle2)), 
                     self.dec + np.amax(self.radius)*np.sin(np.amin(self.angle2))] 
        minRa = min(vertex1[0],vertex2[0],vertex3[0])
        maxRa = max(vertex1[0],vertex2[0],vertex3[0]) 
        minDec = min(vertex1[1],vertex2[1],vertex3[1]) 
        maxDec = max(vertex1[1],vertex2[1],vertex3[1])
        return minRa, maxRa, minDec, maxDec  

#Fake class for linkmap efficiency
class Fake(): 
    def __init__(self, fakeid, objid, mjd, exp, ccd, ra, dec, flux, band):
        self.fakeid = fakeid 
        self.objid = objid 
        self.mjd = mjd  
        self.expnum = exp   
        self.mag = -2.5*np.log10(flux)+31.4 
        self.band = band  
        self.ccd = ccd 
        self.ra = ra  
        self.dec = dec     
        self.links = []
        self.alllinks = []

    def toStr(self): 
        fakeid = "fakeid: " + str(self.fakeid) + "  " 
        objid = "objid: " + str(self.objid) + "  " 
        mjd = "mjd: " + str(self.mjd) + "  "  
        exp = "expnum: " + str(self.expnum) + " "   
        ccd = "ccd: " + str(self.ccd) + " " 
        ra = "RA: " + str(self.ra) + ' '  
        dec = "Dec: " + str(self.dec) + ' '  
        mag = "Mag: " + str(self.mag) + ' ' 
        band = "Band: " + str(self.band)  
        return (fakeid + objid + mjd + exp + ccd + ra + dec + mag + band)
    def toDat(self):
        s = (str(self.mjd+2400000.5) + ' ' + str(ephem.hours(self.ra*np.pi/180)) +
                ' ' + str(ephem.degrees(self.dec*np.pi/180)) + ' 0.1 807')
        return s
 

#class that stores triplets
class Triplet:
    #input is array of detections
    def __init__(self, dets):
        #is a list
        self.dets = dets
        #reduce runtime
        self.elements = 0
        self.errs = 0
        self.orbit = 0
        self.chiSq = -1
    #converts triplet to a dat format
    def toDat(self, outfile):
        s = self.dets[0].toDat()
        for x in range(len(self.dets)-1):
            s = s + '\n'
            s = s + self.dets[x+1].toDat()    
        print('writing to file: ' + outfile)
        with open(outfile, 'w+') as f:
            f.write(s)
        return s

    def calcOrbit(self):
        #converts every coordinate into correct units
        #time0 = time.time()
        mag = max([x.mag for x in self.dets])
        errs=0.1
        if mag <= 21:
            errs = 0.1
        else:
            errs = 0.1 + (mag-21.0) / 10.0
        ralist = [ephem.hours(np.deg2rad(det.ra)) for det in self.dets]
        #time1 = time.time()
        declist = [ephem.degrees(np.deg2rad(det.dec)) for det in self.dets] 
        #time2 = time.time()
        datelist = [ephem.date((Time(det.mjd,format='mjd')).datetime) for det in self.dets]
        #time3 = time.time()
        orbit = Orbit(dates=datelist, ra=ralist, dec=declist, 
                obscode=np.ones(len(self.dets), dtype=int)*807, err=errs)
        self.orbit = orbit
        #time4 = time.time()
        #print('time1: ' + str(time1-time0))
        #print('time2: ' + str(time2-time1))
       # print('time3: ' + str(time3-time2))
        #print('time4: ' + str(time4-time3))
        orbit.get_elements()
        self.chiSq = orbit.chisq
        self.elements, self.errs = orbit.get_elements()
        return self.elements, self.errs

    def setOrbit(self):
        mag = max([x.mag for x in self.dets])
        errs=0.1
        if mag <= 21:
            errs = 0.1
        else:
            errs = 0.1 + (mag -21.0) / 10.0
        ralist = [ephem.hours(np.deg2rad(det.ra)) for det in self.dets]
        declist = [ephem.degrees(np.deg2rad(det.dec)) for det in self.dets]
        datelist = [ephem.date((Time(det.mjd,format='mjd')).datetime) for det in self.dets]
        orbit = Orbit(dates=datelist, ra=ralist, dec=declist,
            obscode=np.ones(len(self.dets), dtype=int)*807, err=errs)
        self.orbit = orbit
        return orbit

    def getDistance(self):
        if self.orbit == 0:
            self.setOrbit()
        dist, err = self.orbit.barycentric_distance()
        return dist, err

    #get the covarance matrix 
    def getCovar(self):
        if(self.orbit==0):
            self.setOrbit()
        #_, _, covar = self.orbit.get_elements_abg()
        return self.orbit.covar_abg

    #gets the chi squared of the orbit
    def getChiSq(self):
        if(self.chiSq == -1):
            self.setOrbit()
            self.chiSq = self.orbit.chisq 
        #return self.orbit.chisq
        return self.chiSq
    #predicts the position of the triplet at a future date
    '''
    input: date in the format MJD
    output: a tuple containing (ra, dec) and a scalar error in arc seconds
    '''
    def predictPos(self, date):
        date = ephem.date((Time(date, format='mjd')).datetime)
        if(self.orbit==0):
            self.setOrbit()
        orbit = self.orbit
        predPos = orbit.predict_pos(date)
        ra = ephem.hours(predPos['ra'])
        dec = ephem.degrees(predPos['dec'])
        err = predPos['err']['a']
        ra = np.degrees(ra)
        dec = np.degrees(dec)
        return (ra, dec), err
    #checks if the objects all have the same fakeid    
    def sameFake(triplet):
        id = triplet.dets[0].fakeid
        if(id == 0):
            return False
        for det in triplet.dets: 
            if(det.fakeid != id):
                return False  
        return True 
 
    def sortByMjd(self):
        self.dets.sort(key = lambda x: x.mjd)

    def apartByXDays(triplet, days):
        mjd = sorted([x.mjd for x in triplet.dets])
        for i in range(1, len(mjd)):
            if mjd[i] - mjd[i - 1] < days:
                return False
        return True
    
    #checks if a triplet is small enough to be a p9 candidate
    def p9Like(self):
        apartDays = 7
        apartDist = 90 # in arcseconds
        if(self.apartByXDays(apartDays)):
            self.sortByMjd()            
            for x in range(len(self.dets)-1):
                dist = calcDist(self.dets[x], self.dets[x+1])
                if(dist*60*60 > apartDist):
                    return False

            return True
        else:
            return False
    def __hash__(self):
        return hash(tuple([x.objid for x in self.dets]))

    def __eq__(self, other):
        selfIDs = tuple([x.objid for x in self.dets]) 
        otherIDs = tuple([x.objid for x in other.dets]) 
        return sorted(selfIDs) == sorted(otherIDs)

    def __ne__(self, other):
        selfIDs = tuple([x.objid for x in self.dets])
        otherIDs = tuple([x.objid for x in other.dets])
        return not (sorted(selfIDs) == sorted(otherIDs))

    def __str__(self):
        return self.toStr()
     
    def toStr(self):
        s = '\n *****Triplet*****\n'
        for x in range(len(self.dets)):
            s += ('Det' + str(x) + ': ' + self.dets[x].toStr() + '\n')
        s+= 'elements: ' + str(self.elements) + '\n'
        s+= 'errs: ' + str(self.errs) + '\n'
        return s

#a class for processing detections in the csvfile, organizing them into fakeObjs
class fakeObj:
    def __init__(self, fakeid):
        self.fakeid = fakeid
        self.listobj = []
        self.campaigns = []

#Write the triplet list to a file
def writeTriplets(tripList, outfile, writeOrbit):
    print('\nwriting triplets to: ' + outfile)
    with open(outfile, 'w+') as output:
        count = 0
        time0 = time.time()
        for triplet in tripList:
            #printPercentage(count, len(tripList), time.time() - time0)
            if(triplet.elements == 0 and triplet.errs == 0 and writeOrbit == True):
                elements, errs = triplet.calcOrbit()
                triplet.getChiSq() 
            output.write('triplet ' + str(count) + ': ' + triplet.toStr())
            elements, errs = triplet.elements, triplet.errs
            output.write('chisq: ' + str(triplet.chiSq) + '\n')
            output.write('*****************\n\n')
            count += 1
    print('\ndone writing')

# sets all orbit to 0 before pickling because Orbit cannot be pickled
def pickleTriplets(tripList, outfile):
    for trip in tripList:
        trip.orbit = 0
    with open(outfile, 'wb') as f:
        print('\npickling to ' + outfile)
        pickle.dump(tripList, f)

def approx(a, b, dist):
    return a-b<=dist and b-a <= dist
def printPercentage(progress, finish, time):
    percent = float(progress)/float(finish) *100.0
    sys.stdout.write("\rCompleted: %1f %% after %5f sec" % (percent, time))
    sys.stdout.flush()

#calculates the distance between two detections in degreees
def calcDist(det1, det2):
    deltaRa = det1.ra - det2.ra
    deltaDec = det1.dec - det2.dec
    dist = deltaRa*deltaRa+deltaDec*deltaDec
    return np.sqrt(dist)
#takes a variable and tells whether it's a float
def isNumber(s):
    try: 
        float(s)
        return True
    except ValueError:
        return False
