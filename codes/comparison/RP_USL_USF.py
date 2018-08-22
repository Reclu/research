#!/usr/bin/pyton

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import pdb
import os

directory=os.path.basename(__file__)[:-3]

if not os.path.exists('texFiles/'+str(directory)):
    os.system('mkdir texFiles/'+str(directory))
path='texFiles/'+str(directory)


"""
Comparison of the implementation of the 1D elastic set of equations
in dynamics:
- with the USL
- with the GIMP
"""

def export2DTeXFile(fileName,xFields,xlabel,ylabel,time,yfields,*kwargs):
    TeXFile=open(fileName,"w")
    n_fields = np.shape(yfields)[0]
    n_labels = np.shape(kwargs)[0]
    # Define Paul Tol's colors (purple to red)
    color=['Purple','Orange','Blue','Red','Duck','Green','Yellow']
    col=['120,28,129','231,133,50','63,96,174','217,33,32','83,158,182','109,179,136','202,184,67']
    marker=['*','x','triangle*','square*','+','star','pentagone*']
    marker=['none','none','+','triangle','none','star','pentagone*']
    style=['solid','dashed','solid','solid','solid','star','pentagone*']
    thickness=['very thick','very thick','thin','thin','thick']
    couleur=['Purple','Orange','Blue','Red','black','Duck','Green']
    for i in range(len(col)):
        TeXFile.write(r'\definecolor{'+color[i]+'}{RGB}{'+col[i]+'}')
        TeXFile.write('\n')
    TeXFile.write(r'\begin{tikzpicture}[scale=0.9]');TeXFile.write('\n')
    TeXFile.write(r'\begin{axis}[xlabel='+str(xlabel)+',ylabel='+str(ylabel)+',ymajorgrids=true,xmajorgrids=true,title={() '+time+'}]');TeXFile.write('\n')
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
    TeXFile.write(r'\legend{'+str(legend)+'}')
    TeXFile.write('\n')    
    TeXFile.write(r'\end{axis}')
    TeXFile.write('\n')
    TeXFile.write('\end{tikzpicture}')
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
CFL=0.7
NTmaxi = 300
length = 6.0
ppc=1
Nelem = 50
E = 2.0e11
Sigy = 400.0e7
H = 10e9
rho = 7800.0
c=np.sqrt(E/rho)
sigd =0.# -0.25*Sigy
v0=0.5*Sigy/(2*rho*c)
factor=1.
timeOut = 0.36*length/(np.sqrt(E/rho))#0.002
t_order=1
timeUnload = 2*timeOut#2.e-4#2*timeOut
algo = 'USL'
# limit = 0 : minmod // limit = 1 : superbee // limit = 2 : MUSCL
limit=-1
update_position=False
mpm_mapping=True

parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy, "H":H,"rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"limit":limit,"algo":algo,"t_order":t_order,"mpm_mapping":mpm_mapping}
#################


##MPM: Material Point Method
USL = dict(parameters)
print 'Computing MPM'
execfile('mpm/elasticity.py', USL)

algo='USF'
parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy, "H":H,"rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"limit":limit,"algo":algo,"t_order":t_order,"mpm_mapping":mpm_mapping}
##MPM: Material Point Method
USF = dict(parameters)
print 'Computing modified MPM'
execfile('mpm/elasticity.py', USF)


algo='USL'
ppc=2
parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy, "H":H,"rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"limit":limit,"algo":algo,"t_order":t_order,"mpm_mapping":mpm_mapping}

print "=============== 2PPC COMPUTATIONS ===================="

##MPM: Material Point Method
USL2 = dict(parameters)
print 'Computing USL'
execfile('mpm/elasticity.py', USL2)

algo='USF'
parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy, "H":H,"rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"limit":limit,"algo":algo,"t_order":t_order,"mpm_mapping":mpm_mapping}
##MPM: Material Point Method
USF2 = dict(parameters)
print 'Computing modified MPM'
execfile('mpm/elasticity.py', USF2)


#############################################################################
#########################  Comparison  ######################################
#############################################################################
x = np.linspace(0.,length,Nelem+1)
dx = x[1]-x[0]
y=x+(dx/2.0)
#Centres of cells/Elements
y=y[:(len(y)-1)]

####Animated plot ###########################################################
from matplotlib import rcParams
rcParams['axes.labelsize'] = 16
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] = 16
rcParams['legend.fontsize'] = 16

#frames=[20,30,40,50,60,80,100,120,140]
frames=[10,25]
#pdb.set_trace()
c = np.sqrt(E/rho)
HT = (H*E)/(H+E) ; cp = np.sqrt(HT/rho)

## Middle point velocity plots

plt.plot(USL["time"][:-1],USL["velo"][24,:],label="USL 1ppc")
plt.plot(USL2["time"][:-1],USL2["velo"][49,:],label="USL 2ppc")
plt.plot(USF["time"][:-1],USF["velo"][24,:],label="USF 1ppc")
plt.plot(USF2["time"][:-1],USF2["velo"][49,:],label="USF 2ppc")
plt.grid()
plt.legend()
plt.show()

# plt.plot(USF["time"][:-1],USF["NRG"][:-1]/max(USF["NRG"][:-1]))
# plt.grid()
# plt.show()

plt.plot(USL["time"][:-1],USL["NRG"][:-1]/max(USL["NRG"][:-1]),'b-x',lw=2.,label='USL 1ppc')
plt.plot(USL2["time"][:-1],USL2["NRG"][:-1]/max(USL2["NRG"][:-1]),'r-x',lw=2.,label='USL 2ppc')
plt.plot(USF["time"][:-1],USF["NRG"][:-1]/max(USF["NRG"][:-1]),'bo',lw=2.,label='USF 1ppc')
plt.plot(USF2["time"][:-1],USF2["NRG"][:-1]/max(USF2["NRG"][:-1]),'ro',lw=2.,label='USF 2ppc')
plt.grid()
plt.legend(numpoints=1)
plt.show()

export2DTeXFile(str(path)+'/US_energies.tex',np.array([USL["time"][:-1],USL2["time"][:-1],USF["time"][:-1],USF2["time"][:-1]]),'$time (s)$',"$\frac{e}{e_{max}}$",'Evolution of total energy',np.array([USL["NRG"][:-1]/max(USL["NRG"][:-1]),USL2["NRG"][:-1]/max(USL2["NRG"][:-1]),USF["NRG"][:-1]/max(USF["NRG"][:-1]),USF2["NRG"][:-1]/max(USF2["NRG"][:-1])]),['USL 1ppc','USL 2ppc','USF 1ppc','USF 2ppc'])

N2=int(Nelem/2.)-40
N1=int(Nelem/2.)-50

for n1 in frames:
    time = '%.2e' % USL["time"][n1]
    plt.plot(USL["pos"][:,n1],USL["Sth"][:,n1],'k-',lw=2.,ms=8.,label='analytical')
    plt.plot(USL["pos"][:,n1],USL["sig"][:,n1],'g-x',lw=2.,ms=8.,label='USL 1ppc')
    plt.plot(USL2["pos"][:,n1],USL2["sig"][:,n1],'g-o',lw=2.,ms=8.,label='USL 2ppc')
    plt.plot(USF["pos"][:,n1],USF["sig"][:,n1],'rx',lw=2.,ms=8.,label='USF 1ppc')
    plt.plot(USF2["pos"][:,n1],USF2["sig"][:,n1],'ro',lw=2.,ms=8.,label='USF 2ppc')
    plt.title('Contrainte longitudinale dans la barre au temps t='+str(time)+' s.',size=24.)
    plt.xlabel('x (m)',size=24.)
    plt.ylabel(r'$\sigma (Pa)$',size=28.)
    plt.legend(numpoints=1)
    plt.grid()
    plt.show()
    
    export2DTeXFile(str(path)+'/US_diffusion'+str(n1)+'.tex',np.array([USL["pos"][:,n1],USL2["pos"][:,n1],USF["pos"][:,n1],USF2["pos"][:,n1]]),'$x (m)$','$\sigma (Pa)$',str(time),np.array([USL["sig"][:,n1],USL2["sig"][:,n1],USF["sig"][:,n1],USF2["sig"][:,n1]]),['USL 1ppc','USL 2ppc','USF 1ppc','USF 2ppc'])
    export2DTeXFile(str(path)+'/US_velo'+str(n1)+'.tex',np.array([USL["pos"][:,n1],USL2["pos"][:,n1],USF["pos"][:,n1],USF2["pos"][:,n1],USF2["pos"][:,n1]]),'$x (m)$','$v (m/s)$',str(time),np.array([USL["velo"][:,n1],USL2["velo"][:,n1],USF["velo"][:,n1],USF2["velo"][:,n1],USF2["Vth"][:,n1]]),['USL 1ppc','USL 2ppc','USF 1ppc','USF 2ppc','analytical'])
