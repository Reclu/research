#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import pdb

"""
Comparison between two implementation of DGMPM.
One aims to reduce numerical diffusion by interpolating variation from nodes two material points after resolution of the discrete equations on the mesh
"""

ppc=4
Mp=100*ppc
if ppc==1: CFL = 0.5
elif ppc==2: CFL=0.4286
elif ppc==4: CFL=0.2258
parameters={"algo":'test',"sinusoidal":False,"Mp":Mp,"ppc":ppc,"CFL":CFL}
DGMPM = dict(parameters)
execfile('dgmpm.py', DGMPM)

parameters={"algo":'original',"sinusoidal":False,"Mp":Mp,"ppc":ppc,"CFL":CFL}
DGMPM2 = dict(parameters)
execfile('dgmpm.py', DGMPM2)

DGFEM = dict(parameters)
execfile('dgfem.py', DGFEM)

fig = plt.figure()

plt.grid()
plt.xlim(0.,1.)
plt.ylim(-100.,100.)
plt.xlabel('x (m)', fontsize=18)
plt.ylabel(r'$\mathcal{Q}$ (Pa)', fontsize=18)
plt.title('Stress wave propagation in a bar', fontsize=16)
line1,= plt.plot([], [],'r-o', lw=1.5)
line2,= plt.plot([], [],'b-s', lw=1.5)
line3,= plt.plot([], [],'k', lw=1.5)
fig.legend((line1,line2,line3),('DGMPM test','DGMPM original','Analytical'),'upper right',numpoints=1)

# initialization function: plot the background of each frame
def init():
    line2.set_data([], [])
    line3.set_data([], [])
    return line2,line3

# animation function.  This is called sequentially
def animate(i):
    line1.set_data(DGMPM["xp"][:,0],DGMPM["Stress"][:,i])
    line2.set_data(DGMPM2["xp"][:,0],DGMPM2["Stress"][:,i])
    line3.set_data(DGMPM2["xp"][:,0],DGMPM["analytical"][:,i])
    return line1,line2,line3
 
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=DGMPM2["Stress"].shape[1], interval=50, blit=True)

#Animation of the stress
plt.show()
#anim.save('StressBar.mp4', extra_args=['-vcodec', 'libx264'])

