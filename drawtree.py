'''
Draw a merger tree
Sam Geen, December 2016
'''

import numpy as np
import matplotlib.pyplot as plt
import readclumps, trackclumps

def Draw(tree):
    tracks = tree.Tracks()
    for track in tracks.itervalues():
        pos = np.array(track.Position())
        t = track.Time()
        x = pos[:,0]
        y = pos[:,1]
        z = pos[:,2]
        plt.plot(x,y)
    plt.savefig("plots/test.pdf")

if __name__=="__main__":
    loc = "/data/Simulations/Clumps/run2703/clumps"
    sim = readclumps.Simulation(loc)
    tree = trackclumps.ClumpTree(sim)
    Draw(tree)
    
