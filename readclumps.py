'''
Reads clumps from a timestep
Sam Geen, December 2016
'''

import glob, os
import numpy as np

from collections import OrderedDict

class Snapshot(object):
    '''
    A snapshot of the simulation with all clumps and sinks
    TODO: Split off sinks into new object
    '''
    def __init__(self,sim,step):
        self._sim = sim
        self._step = step
        self._info = None
        self._clumps = None
        self._clumpheader = None
        self._props = None
        self._deriveds = ["time","radius","pos","vel"]
        # Set up
        # TODO: Do this on-the-fly to make faster
        self._ReadInfo()
        self._Read()

    def __getitem__(self,prop):
        '''
        Get all the properties with this value
        '''
        if prop in self._deriveds:
            return self._Derived(prop)
        return self._clumps[:,self._props[prop]]

    def _Derived(self,prop):
        '''
        Return a derived value
        '''
        if prop == "time":
            # Time of main step
            return self._info["time"]
        if prop == "radius":
            # R = (3 M / 4 pi rho_av)^(1/3)
            return (3.0 / (4.0 * np.pi) * 
                    self["mass_cl"] / self["rho_av"])**(1.0/3.0)
        if prop == "pos":
            # Concatenate position
            return np.array([self["peak_x"],self["peak_y"],self["peak_z"]])
        if prop == "vel":
            # Concatenate velocity
            return np.array([self["v_x"],self["v_y"],self["v_z"]])
        print "No derived property with the name", prop
        raise ValueError

    def _Read(self):
        search = self._sim.Folder()+"/clump_"+str(self._step).zfill(5)+"*"
        files = glob.glob(search)
        readme = self._sim.Folder()+"/README"
        self._clumpheader = open(readme).readlines()[1].split()
        self._props = {p:i for p, i in zip(self._clumpheader,
                                           range(0,len(self._clumpheader)))}
        data = None
        for f in files:
            d = np.atleast_2d(np.loadtxt(f))
            if data is None:
                data = d
            else:
                try:
                    data = np.concatenate((data,d))
                except:
                    print "ERRAR"
                    import pdb; pdb.set_trace()
        # Sort by x position to make matching easier
        data = data[data[:,self._props["mass_cl"]].argsort()[::-1]]
        self._clumps = data

    def _ReadInfo(self):
        self._info = {}
        # Read info file
        fname = self._sim.Folder()+"/info_"+str(self._step).zfill(5)+".dat"
        lines = open(fname).readlines()
        # Time
        for line in lines:
            lhs, rhs = line.split("=")
            self._info[lhs.strip()] = float(rhs.strip())

class Simulation(object):
    def __init__(self,folder):
        self._folder = folder
        self._allsteps = None
    
    def Folder(self):
        return self._folder

    def FindSnapshot(self,step):
        '''
        Returns a snapshot with the given step number
        '''
        step = int(step)
        if not step in self.AllSteps():
            print "Step",step,"not in simulation!"
            raise ValueError
        return Snapshot(self,step)
        
    def AllSteps(self):
        '''
        Get a list of all step numbers
        '''
        # Request steps if not set yet
        if self._allsteps is None:
            # Search for any sink/clump files
            search = self._folder+"/clump_*"
            files = glob.glob(search)
            self._folder+"/sink_*"
            files += glob.glob(search)
            # Get list of unique step numbers in order
            nums = [int(file[len(search)-1:len(search)+4]) for file in files]
            nums = set(nums)
            nums = list(nums)
            nums.sort()
            self._allsteps = np.array(nums)
        return self._allsteps



if __name__=="__main__":
    sim = Simulation("/data/Simulations/Clumps/run2703/clumps")
    snaps = []
    for step in sim.AllSteps():
        snaps.append(sim.FindSnapshot(step))
    print "DUN"
