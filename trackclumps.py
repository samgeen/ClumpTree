'''
Tracks clumps between timesteps
Sam Geen, December 2016
'''
import numpy as np

import readclumps

def mag(vect):
    '''
    Find the magnitude of array of vectors
    arr(3,N) -> arr(N)
    '''
    return np.sqrt(np.sum(vect*vect,axis=0))

class Track(object):
    '''
    Stores a track for a given clump over time
    '''
    # NOTE: id is a python keyword
    def __init__(self,myid):
        self._id = myid
        self._time = []
        self._pos = []
        self._mass = []
        self._radius = []
        
    def Add(self,time,pos,mass,radius):
        self._time.append(time)
        self._pos.append(pos)
        self._mass.append(mass)
        self._radius.append(radius)

    def Time(self):
        return self._time

    def Position(self):
        return self._pos

    def Mass(self):
        return self._mass

    def Radius(self):
        return self._radius

class ClumpTree(object):
    '''
    Merger tree for clumps
    '''
    def __init__(self,sim):
        self._sim = sim
        self._snaps = None
        # Clump data tracks
        self._tracks = None
        # Matches one timestep to the next one
        self._matcher = Matcher()
        # TODO: REFACTOR THIS OUT OF THE CONSTRUCTOR
        self._MakeTree()
        # Current unique ID
        self._maxid = self._matcher.MaxID()

    def Sim(self):
        return self._sim

    def IDs(self):
        return self._snaps

    def Tracks(self):
        return self._tracks

    def _MakeTree(self):
        # Find the snapshot IDs for each step
        steps = self._sim.AllSteps()
        # Fill first snapshot's clumps
        clumps = self._sim.FindSnapshot(steps[0])
        snaps = [self._matcher.Start(clumps)]
        # Match each successive snapshot
        tracks = {}
        for step in steps[1:]:
            # Find IDs
            clumps = self._sim.FindSnapshot(step)
            ids = self._matcher.FindNextIDs(clumps)
            # Make tracks for each clump ID
            for ic,cid in zip(range(0,len(ids)),ids):
                # Add current ID if not in track
                if cid not in tracks:
                    tracks[cid] = Track(cid)
                time = clumps["time"]
                pos = clumps["pos"][:,ic]
                mass = clumps["mass_cl"][ic]
                radius = clumps["radius"][ic]
                tracks[cid].Add(time,pos,mass,radius)
            # Populate clump tracks
            snaps.append(ids)
        self._snaps = snaps
        self._tracks = tracks

class Matcher(object):
    '''
    Class that matches two snapshots and passes out IDs
    '''
    def __int__(self):
        # Old and new clumps
        self._cuid = 0
        self._cold = None
        self._cnew = None
        self._ids = None

    def MaxID(self):
        '''
        Largest particle ID
        '''
        return self._cuid
        
    def Start(self,cfirst):
        '''
        Populate first step with IDs
        '''
        lf = len(cfirst["index"])
        self._ids = np.arange(0,lf)+1
        self._cuid = lf+1
        self._cold = cfirst
        return self._ids

    def FindNextIDs(self,cnew):
        '''
        Find next IDs from new set of clumps
        '''
        # More efficient to re-use class???
        self._cnew = cnew
        # Find IDs from one step to the next
        oldids = self._ids
        # Matrix containing quality of fit with new clumps
        quality = self._MatchClumps()
        # Assign IDs based on quality matrix
        ids = self._AssignIDs(quality)
        # Clean up
        self._cold = cnew
        self._ids = ids
        # Return
        return self._ids

    def _MatchClumps(self):
        # Find a matrix of quality-of-fit between clumps
        # Each value should be 0 (terrible fit) to 1 (identical)
        # This is a combination of values
        # Each value is 1 - a*f, where 
        # f is the fractional difference between two quantities, and
        # a is some scaling factor
        # Set up comparison
        c0 = self._cold
        c1 = self._cnew
        dt = c1["time"] - c0["time"]
        n0 = len(c0["index"])
        n1 = len(c1["index"])
        # Some convenience definitions
        p0 = c0["pos"]
        p1 = c1["pos"]
        r0 = c0["radius"]
        r1 = c1["radius"]
        m0 = c0["mass_cl"]
        m1 = c1["mass_cl"]
        # DRIFT DISTANCE
        # Compare drifted position p1 - (p0 + v0)
        pdrift = c0["pos"]+c0["vel"]*dt
        # TODO: optimise for speed somehow
        f = np.zeros((n0,n1))
        # This is the meat of the matching algorithm
        # Heuristic as all hell, beware
        for i0 in range(0,n0):
            for i1 in range(0,n1):
                # Heuristic distance measure
                # Allow up to 3*r1 drift before failing fit
                ddrift = mag(p1[:,i1]-pdrift[:,i0])
                fdrift = max(0,1.0 - 0.3*ddrift / r1[i1])
                # Radius (allow 3x change before failing)
                dradius = np.abs(r1[i1] - r0[i0]) / (0.5 * (r1[i1] + r0[i0]))
                fradius = max(0,1.0 - 0.3*dradius)
                # Mass (allow 3x change before failing)
                dmass = np.abs(m1[i1] - m0[i0]) / (0.5 * (m1[i1] + m0[i0]))
                fmass = max(0,1.0 - 0.3*dradius)
                # Combine heuristics
                f[i0,i1] = fdrift*fradius*fmass
        return f

    def _NewID(self):
        '''
        Get a new ID from the shelf
        '''
        # Increment current unique ID
        self._cuid += 1
        return self._cuid

    def _AssignIDs(self,quality):
        '''
        Assign IDs based on the quality of fit matrix
        '''
        ids0 = self._ids
        ids1 = []
        # For each new clump, figure out best old clump
        # If no match, make a new ID
        n0, n1 = quality.shape
        # New clump should be sorted in descending mass order
        # Every time we assign an old clump to a new one,
        # empty out the old IDs
        for i1 in range(0,n1):
            q = quality[:,i1]
            ind = np.argmax(q)
            best = q[ind]
            newid = ids0[ind]
            # If no match, assign new ID (new clump!)
            if best == 0:
                newid = self._NewID()
            # Unassign chosen ID to prevent multiple assignation
            # TODO: Allow sub-clumps
            q[ind] = 0.0
            ids1.append(newid)
        return ids1

if __name__=="__main__":
    loc = "/data/Simulations/Clumps/run2703/clumps"
    sim = readclumps.Simulation(loc)
    tree = ClumpTree(sim)
    print tree.IDs()
