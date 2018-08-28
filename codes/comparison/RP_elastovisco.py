#!/usr/bin/pyton

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import pdb
import os
import sys
directory=os.path.basename(__file__)[:-3]

if not os.path.exists('texFiles/'+str(directory)):
    os.system('mkdir texFiles/'+str(directory))
path='texFiles/'+str(directory)
"""
Comparison of the implementation of the 1D elastic set of equations
in dynamics:
- with the MPM
- with the GIMP
"""
def export2DTeXFile(fileName,xFields,xlabel,ylabel,subtitle,yfields,*kwargs):
    TeXFile=open(fileName,"w")
    n_fields = np.shape(yfields)[0]
    n_labels = np.shape(kwargs)[0]
    # Define Paul Tol's colors (purple to red)
    marker=['+','x','triangle','square','none','none','pentagone*']
    #marker=['none','none','+','triangle','none','star','pentagone*']
    style=['only marks','only marks','only marks','only marks','dashed','solid','pentagone*']
    thickness=['very thick','very thick','very thick','very thick','very thick','thin','thin','thick']
    couleur=['Red','Orange','Blue','Purple','Green','black','Duck','Green']
    TeXFile.write(r'\begin{tikzpicture}[scale=0.9]');TeXFile.write('\n')
    TeXFile.write(r'\begin{axis}[xlabel='+str(xlabel)+',ylabel='+str(ylabel)+',ymajorgrids=true,xmajorgrids=true,legend pos=outer north east,title={'+subtitle+'}]');TeXFile.write('\n')
    legend=''
    for i in range(n_fields):
        if i==0:
            legend=legend+kwargs[0][i]
        else:
            legend=legend+','+kwargs[0][i]
        TeXFile.write(r'\addplot['+str(couleur[i])+','+str(thickness[i])+',mark='+str(marker[i])+','+str(style[i])+'] coordinates {')
        for j in range(np.shape(yfields[i])[0]):
            TeXFile.write('('+str(xFields[i][j])+','+str(yfields[i][j])+') ')
        TeXFile.write('};\n')
    if subtitle=='(c) evolution of total energy $e$':
        TeXFile.write(r'\legend{'+str(legend)+'}')
    else:
        TeXFile.write(r'%\legend{'+str(legend)+'}')
    TeXFile.write('\n')    
    TeXFile.write(r'\end{axis}')
    TeXFile.write('\n')
    TeXFile.write('\end{tikzpicture}')
    TeXFile.write('\n')
    TeXFile.write('%%% Local Variables:')
    TeXFile.write('\n')
    TeXFile.write('%%% mode: latex')
    TeXFile.write('\n')
    TeXFile.write('%%% TeX-master: "../../mainManuscript"')
    TeXFile.write('\n')
    TeXFile.write('%%% End:')
    TeXFile.write('\n')
    TeXFile.close()


def export2pgfPlot(fileName,xfield,yfield,xlabel,ylabel):
    #pdb.set_trace()
    dataFile=open(fileName,"w")
    dataFile.write('# Curve ('+str(xlabel)+';'+str(ylabel)+') '+str(len(xfield))+' points.\n')
    for i,x in enumerate(xfield):
        dataFile.write(str(x)+' '+str(yfield[i])+' i\n')
    dataFile.close()

###Opening the files and computation of the solution by each method
###Parameters####
CFL=0.5
NTmaxi = 300
length = 6.0
ppc=1
Nelem = 50
E = 2.0e11
Sigy = 400.0e6
H = 10e9
rho = 7800.0
eta=100.0
n=1./4.
c=np.sqrt(E/rho)
sigd =0.
v0=5.*Sigy/(rho*c)
hardening='isotropic'
factor=1.
timeOut = 1.*length/(np.sqrt(E/rho))
t_order=1
timeUnload = 2*timeOut
algo = 'USL'
update_position=False
mpm_mapping=True
limit=-1

parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy, "H":H,"rho":rho,"eta":eta,"n":n,"sigd":sigd,"hardening":hardening,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"algo":algo,"t_order":t_order,"limit":limit,"mpm_mapping":mpm_mapping,"compute_CFL":False}
#################


##MPM: Material Point Method
MPM = dict(parameters)
print 'Computing MPM'
execfile('mpm/elastoviscoplasticity.py', MPM)


##DGMPM: Discontinuous Galerkin Material Point Method
DGMPM = dict(parameters)
print 'Computing  DGMPM'
if hardening=='kinematic':
    execfile('dgmpm/EVP_kin.py', DGMPM)
elif hardening=='isotropic':
    execfile('dgmpm/EVP_iso.py', DGMPM)


##FVM: Finite Volume Method
FVM = dict(parameters)
print 'Computing FVM'
execfile('fvm/elastoviscoplasticity.py', FVM)

ppc=2
#Nelem = 50/ppc
parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy, "H":H,"rho":rho,"eta":eta,"n":n,"sigd":sigd,"hardening":hardening,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"algo":algo,"t_order":t_order,"limit":limit,"mpm_mapping":mpm_mapping,"compute_CFL":False}

print "=============== 2PPC COMPUTATIONS ===================="

##MPM: Material Point Method
MPM2 = dict(parameters)
print 'Computing MPM'
execfile('mpm/elastoviscoplasticity.py', MPM2)

##DGMPM: Discontinous Galerkin Material Point Method
DGMPM2 = dict(parameters)
print 'Computing DGMPM'
if hardening=='kinematic':
    execfile('dgmpm/EVP_kin.py', DGMPM2)
elif hardening=='isotropic':
    execfile('dgmpm/EVP_iso.py', DGMPM2)

parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy, "H":H,"rho":rho,"eta":eta,"n":n,"sigd":sigd,"hardening":hardening,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"algo":algo,"t_order":2,"limit":limit,"mpm_mapping":mpm_mapping,"compute_CFL":False}
##DGMPM: Discontinous Galerkin Material Point Method
DGMPM3 = dict(parameters)
print 'Computing DGMPM (RK2)'
if hardening=='kinematic':
    execfile('dgmpm/EVP_kin.py', DGMPM3)
elif hardening=='isotropic':
    execfile('dgmpm/EVP_iso.py', DGMPM3)
#############################################################################
#########################  Comparison  ######################################
#############################################################################

####Animated plot ###########################################################
from matplotlib import rcParams
rcParams['axes.labelsize'] = 16
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] = 16
rcParams['legend.fontsize'] = 16

frames=[1,2,3,4,5,6]
frames=[]
for i in frames:
    plt.plot(MPM["pos"][:,i],MPM["sig"][:,i],'k')
    plt.plot(DGMPM["pos"][:,i],DGMPM["sig"][:,i],'r--')
    plt.grid()
    plt.show()


c = np.sqrt(E/rho)
HT = (H*E)/(H+E) ; cp = np.sqrt(HT/rho)

"""
plt.plot(DGMPM["xp"][:,0],np.zeros(len(DGMPM["xp"][:,0])),'bo')
plt.plot(DGMPM2["xp"][:,0],np.zeros(len(DGMPM2["xp"][:,0])),'ro')
plt.grid()
plt.show()
"""

factor_dgmpms=float(DGMPM2["increments"]/DGMPM["increments"])
factor_mpms=MPM2["increments"]/MPM["increments"]

factor2=DGMPM2["CFL"]/CFL


### Energy plots
"""
plt.plot(MPM["time"][:-1],MPM["NRG"][:-1]/max(MPM["NRG"][:-1]),'b-x',lw=2.,label='mpm 1ppc')
plt.plot(MPM2["time"][:-1],MPM2["NRG"][:-1]/max(MPM2["NRG"][:-1]),'r-x',lw=2.,label='mpm 2ppc')
plt.plot(DGMPM["time"][:-1],DGMPM["NRG"][:-1]/max(DGMPM["NRG"][:-1]),'bo',lw=2.,label='dgmpm 1ppc')
plt.plot(DGMPM2["time"][:-1],DGMPM2["NRG"][:-1]/max(DGMPM2["NRG"][:-1]),'ro',lw=2.,label='dgmpm 2ppc')
plt.grid()
plt.legend(numpoints=1)
plt.show()

legend=['mpm 1ppc','mpm 2ppc','dgmpm 1ppc','dgmpm 2ppc','dgmpm 2ppc (RK2)','exact']
export2DTeXFile(str(path)+'/dgmpm_mpm_energies.tex',np.array([MPM["time"][:-1],MPM2["time"][:-1],DGMPM["time"][:-1],DGMPM2["time"][:-1],DGMPM3["time"][:-1],np.array([0.,1.e-8])]),'$time (s)$',r'$\frac{e}{e_{max}}$','(c) evolution of total energy $e$',np.array([MPM["NRG"][:-1]/max(MPM["NRG"][:-1]),MPM2["NRG"][:-1]/max(MPM2["NRG"][:-1]),DGMPM["NRG"][:-1]/max(DGMPM["NRG"][:-1]),DGMPM2["NRG"][:-1]/max(DGMPM2["NRG"][:-1]),DGMPM3["NRG"][:-1]/max(DGMPM3["NRG"][:-1]),np.array([1.,1.])]),legend)
"""

frames=[5,20]
frames=[]
for n1 in frames:
    time = '%.2e' % MPM["time"][2*n1]
    plt.plot(MPM["pos"][:,n1],MPM["sig"][:,2*n1],'g-x',lw=2.,ms=8.,label='MPM 1ppc')
    plt.plot(MPM2["pos"][:,n1],MPM2["sig"][:,2*n1],'g-o',lw=2.,ms=8.,label='MPM 2ppc')
    plt.plot(DGMPM["pos"][:,n1],DGMPM["sig"][:,n1],'rx',lw=2.,ms=8.,label='DGMPM 1ppc')
    plt.plot(DGMPM2["pos"][:,n1],DGMPM2["sig"][:,2*n1],'ro',lw=2.,ms=8.,label='DGMPM 2ppc')
    plt.plot(DGMPM3["pos"][:,n1],DGMPM3["sig"][:,n1],'yo',lw=2.,ms=8.,label='DGMPM 2ppc (RK2)')
    plt.title('Contrainte longitudinale dans la barre au temps t='+str(time)+' s.',size=24.)
    plt.xlabel('x (m)',size=24.)
    plt.ylabel(r'$\sigma (Pa)$',size=28.)
    plt.legend(numpoints=1)
    plt.grid()
    plt.show()
    legend=['mpm 1ppc','mpm 2ppc','dgmpm 1ppc','dgmpm 2ppc','dgmpm 2ppc (RK2)']
    if n1==5 : subtitle='(a) time t = '+str(time)+' s.'
    if n1==20 : subtitle='(b) time t = '+str(time)+' s.'
    export2DTeXFile(str(path)+'/dgmpm_mpm_diffusion'+str(n1)+'.tex',np.array([MPM["pos"][:,2*n1],MPM2["pos"][:,2*n1],DGMPM["pos"][:,n1],DGMPM2["pos"][:,2*n1],DGMPM3["pos"][:,n1]]),'$x (m)$',r'$\sigma (Pa)$',str(subtitle),np.array([MPM["sig"][:,2*n1],MPM2["sig"][:,2*n1],DGMPM["sig"][:,n1],DGMPM2["sig"][:,2*n1],DGMPM3["sig"][:,n1]]),legend)
    export2DTeXFile(str(path)+'/dgmpm_mpm_velo'+str(n1)+'.tex',np.array([MPM["pos"][:,2*n1],MPM2["pos"][:,2*n1],DGMPM["pos"][:,n1],DGMPM2["pos"][:,2*n1],DGMPM3["pos"][:,n1]]),'$x (m)$','$v (m/s)$',str(subtitle),np.array([MPM["velo"][:,2*n1],MPM2["velo"][:,2*n1],DGMPM["velo"][:,n1],DGMPM2["velo"][:,2*n1],DGMPM3["velo"][:,n1]]),legend)


print np.shape(DGMPM["pos"]),np.shape(DGMPM2["pos"]),np.shape(MPM2["pos"])
####################################################################
fig, (ax1, ax2) = plt.subplots(2,1)

# intialize two line objects (one in each axes)
line1, = ax1.plot([], [],'r-o', ms=2.5,label='dgmpm')
line2, = ax1.plot([], [],'g--' ,ms=1.5,label='mpm')
line3, = ax1.plot([], [],'rs' , ms=2,label='dgmpm')
line4, = ax1.plot([], [],'b--s' , ms=2,label='mpm')
line5, = ax1.plot([], [],'k' , ms=2,label='fvm')
line6, = ax2.plot([], [],'r-o', ms=2.5,label='dgmpm')
line7, = ax2.plot([], [],'g--' ,ms=1.5,label='mpm')
line8, = ax2.plot([], [],'rs' , ms=2,label='dgmpm')
line9, = ax2.plot([], [],'b--s' , ms=2,label='mpm')
line10, = ax2.plot([], [],'k' , ms=2,label='fvm')
line = [line1, line2,line3,line4,line5,line6,line7,line8,line9,line10]

ax1.grid()
ax2.grid()
ax1.set_xlabel('x (m)', fontsize=18)
ax2.set_xlabel('x (m)', fontsize=18)
ax2.set_ylabel(r'$\varepsilon^p$', fontsize=18)
ax1.set_ylabel(r'$\sigma$', fontsize=18)

ax1.legend(numpoints=1)

ax1.set_xlim(0.,length)
ax2.set_xlim(0.,length)
ax1.set_ylim(1.1*np.min(MPM["sig"]),1.1*np.max(MPM["sig"]))
if hardening=='kinematic':
    ax2.set_ylim(1.1*np.min(MPM2["epsp"]),1.1*np.max(MPM2["epsp"]))
elif hardening=='isotropic':
    ax2.set_ylim(1.1*np.min(MPM2["p"]),1.1*np.max(MPM2["p"]))
ax2.set_ylim(1.1*np.min(MPM2["epsp"]),1.1*np.max(MPM2["epsp"]))
def init():
    line[0].set_data([], [])
    line[1].set_data([], [])
    line[2].set_data([], [])
    line[3].set_data([], [])
    line[4].set_data([], [])
    line[5].set_data([], [])
    line[6].set_data([], [])
    line[7].set_data([], [])
    line[8].set_data([], [])
    line[9].set_data([], [])
    return line

def animate(i):
    line[0].set_data(DGMPM["pos"][:,i],DGMPM["sig"][:,i])
    line[1].set_data(MPM["pos"][:,2*i],MPM["sig"][:,2*i])
    line[2].set_data(DGMPM2["pos"][:,2*i],DGMPM2["sig"][:,2*i])
    line[3].set_data(MPM2["pos"][:,2*i],MPM2["sig"][:,2*i])
    line[4].set_data(FVM["centroids"],FVM["sig"][:,i])
    line[5].set_data(DGMPM["pos"][:,i],DGMPM["epsp"][:,i])
    line[6].set_data(MPM["pos"][:,2*i],MPM["epsp"][:,2*i])
    line[7].set_data(DGMPM2["pos"][:,2*i],DGMPM2["epsp"][:,2*i])
    line[8].set_data(MPM2["pos"][:,2*i],MPM2["epsp"][:,2*i])
    line[9].set_data(FVM["centroids"],FVM["EP"][:,i])
    return line

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=DGMPM["increments"], interval=100, blit=True)

plt.show()
