#!/usr/bin/pyton

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import pdb
import os
import sys
from computeAnalyticalSolutionKIN import *
directory=os.path.basename(__file__)[:-3]

if not os.path.exists('texFiles/'+str(directory)):
    os.system('mkdir texFiles/'+str(directory))
path='texFiles/'+str(directory)
"""
Comparison of the implementation of the 1D elastic set of equations
in dynamics:
- with the MPM
- with the DGMPM
"""
def export2DTeXFile(fileName,xFields,xlabel,ylabel,subtitle,yfields,*kwargs):
    TeXFile=open(fileName,"w")
    n_fields = np.shape(yfields)[0]
    n_labels = np.shape(kwargs)[0]
    # Define Paul Tol's colors (purple to red)
    marker=['+','none','none','x','square','none','none','pentagone*']
    marker=['none','none','+','x','none','none','star','pentagone*']
    style=['dashed','solid','solid','solid','solid','dashed','solid','pentagone*']
    thickness=['very thick','very thick','very thick','thick','thin','very thick','very thick','thin','thin','thick']
    couleur=['Red','Blue','Orange','Purple','black','Orange','Green','Duck','Green']
    TeXFile.write(r'\begin{tikzpicture}[scale=0.9]');TeXFile.write('\n')
    
    if fileName[-10:-6]=='rare':
        TeXFile.write(r'\begin{axis}[xlabel='+str(xlabel)+',ylabel='+str(ylabel)+',ymajorgrids=true,xmajorgrids=true,legend pos= north west,title={'+subtitle+'},xmin=0.,xmax=6.]');TeXFile.write('\n')
    else:
        TeXFile.write(r'\begin{axis}[xlabel='+str(xlabel)+',ylabel='+str(ylabel)+',ymajorgrids=true,xmajorgrids=true,legend pos= south west,title={'+subtitle+'},xmin=0.,xmax=6.]');TeXFile.write('\n')
        
    legend=''
    for i in range(n_fields):
        if i==0:
            legend=legend+kwargs[0][i]
        else:
            legend=legend+','+kwargs[0][i]
        TeXFile.write(r'\addplot['+str(couleur[i])+','+str(thickness[i])+',mark='+str(marker[i])+','+str(style[i])+',mark repeat=5] coordinates {')
        for j in range(np.shape(yfields[i])[0]):
            TeXFile.write('('+str(xFields[i][j])+','+str(yfields[i][j])+') ')
        TeXFile.write('};\n')
    if subtitle[:3]=='(b)':
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

def export2DGroupplot(fileName,containers,rowFields,colFields,titles,Ylabels,legend,*kwargs):
    
    row=len(rowFields)
    col=len(colFields)
    fields_in_plots=len(containers)
    marker=['none','none','+','x','none','none','star','pentagone*']
    style=['dashed','solid','solid','solid','solid','dashed','solid','pentagone*']
    thickness=['very thick','very thick','very thick','thick','thin','very thick','very thick','thin','thin','thick']
    couleur=['Red','Blue','Orange','Purple','black','Orange','Green','Duck','Green']
    maxim=[]
    minim=[]
    maximum=np.zeros(row)
    minimum=np.zeros(row)
    for i,field in enumerate(rowFields):
        for j in range(fields_in_plots):
            for k in (colFields[i]):
                maxim.append(max(containers[j][field][:,k]))
                minim.append(min(containers[j][field][:,k]))
        maximum[i]=1.05*max(maxim)
        minimum[i]=1.05*min(minim)
    TeXFile=open(fileName,"w")
    # Define Paul Tol's colors (purple to red)
    TeXFile.write(r'\begin{tikzpicture}[spy using outlines={rectangle, magnification=3, size=1.5cm, connect spies},scale=.9]');TeXFile.write('\n')
    TeXFile.write(r'\begin{groupplot}[group style={group size='+str(col)+' by '+str(row)+',');TeXFile.write('\n')
    TeXFile.write('ylabels at=edge left, yticklabels at=edge left,horizontal sep=2.ex,');TeXFile.write('\n')
    TeXFile.write('vertical sep=4ex,xticklabels at=edge bottom,xlabels at=edge bottom},');TeXFile.write('\n')
    if row==1:
        TeXFile.write(r'ymajorgrids=true,xmajorgrids=true,enlargelimits=0,xmin=0.,xmax=6.,xlabel=$x (m)$,');TeXFile.write('\n')
    else:
        TeXFile.write(r'ymajorgrids=true,xmajorgrids=true,enlargelimits=0,xmin=0.,xmax=6.,xlabel=$x (m)$,');TeXFile.write('\n')
    TeXFile.write('axis on top,scale only axis,width=0.45\linewidth');TeXFile.write('\n')
    TeXFile.write(']');TeXFile.write('\n')
    for i,field in enumerate(rowFields): ## sum over rows
        for j in range(col):
            TeXFile.write(r'\nextgroupplot[')
            if i==0: TeXFile.write(r'title={'+str(titles[j])+'},ymin='+str(minimum[i])+',ymax='+str(maximum[i])+',')
            if j==0: TeXFile.write(r'ylabel='+str(Ylabels[i])+',ymin='+str(minimum[i])+',ymax='+str(maximum[i])+',')
            if j==col-1 and i==row-1: TeXFile.write(r'legend style={at={($(0.5,-0.35)+(0.45cm,1cm)$)},legend columns=3},ymin='+str(minimum[i])+',ymax='+str(maximum[i]))
            TeXFile.write(']');TeXFile.write('\n')
            for k in range(fields_in_plots):
                TeXFile.write(r'\addplot['+str(couleur[k])+','+str(style[k])+',mark='+str(marker[k])+','+thickness[k]+',mark size=3pt,mark repeat=8] coordinates{')
                #pdb.set_trace()
                #print field
                FIELD=containers[k][field][:,colFields[j][k]]
                xFields=containers[k]["pos"][:,colFields[j][k]]
                for l in range(len(FIELD)):
                    TeXFile.write('('+str(xFields[l])+','+str(FIELD[l])+') ')
                TeXFile.write('};\n')
            if j==0 and fileName[-9:]=='shock.tex':
                TeXFile.write(r'\begin{scope}');TeXFile.write('\n')
                TeXFile.write(r'\spy[black,size=3cm] on (2.35,4.2) in node [fill=none] at (5.8,4.);');TeXFile.write('\n')
                TeXFile.write(r'\end{scope}');TeXFile.write('\n')
    for lab in legend:
        TeXFile.write(r'\addlegendentry{'+str(lab)+'}');TeXFile.write('\n')
    TeXFile.write('\n')    
    TeXFile.write(r'\end{groupplot}')
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
Nelem = 100
E = 2.0e11
Sigy = 400.0e6
rho = 7800.0
nu = 0.3
lamb = (E*nu)/((1+nu)*(1.-2.*nu))
mu = E/(2.*(1.+nu))
C=lamb+2.*mu
c=np.sqrt(C/rho)

factor=1.
timeOut = 1.*length/c
t_order=1
timeUnload = 2*timeOut

sigd = -0.5*Sigy
v0=0.*Sigy/(rho*c)
algo = 'USL'
update_position=False
mpm_mapping=True
limit=-1

parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy,"C":C, "rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"algo":algo,"t_order":t_order,"limit":limit,"mpm_mapping":mpm_mapping,"compute_CFL":False}
#################


##MPM: Material Point Method
MPM = dict(parameters)
print 'Computing MPM'
execfile('mpm/hyperelasticity.py', MPM)


##DGMPM: Discontinuous Galerkin Material Point Method
DGMPM = dict(parameters)
print 'Computing  DGMPM'
execfile('dgmpm/hyperelasticity.py', DGMPM)


# for i,t in enumerate(DGMPM["time"]):
#     print i,t
ppc=2
#Nelem = 50/ppc
parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy,"C":C ,"rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"algo":algo,"t_order":t_order,"limit":limit,"mpm_mapping":mpm_mapping,"compute_CFL":False}


print "=============== 2PPC COMPUTATIONS ===================="

##MPM: Material Point Method
MPM2 = dict(parameters)
print 'Computing MPM'
#execfile('mpm/hyperelasticity.py', MPM2)

##DGMPM: Discontinous Galerkin Material Point Method
DGMPM2 = dict(parameters)
print 'Computing DGMPM'
execfile('dgmpm/hyperelasticity.py', DGMPM2)

parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy,"C":C ,"rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"algo":algo,"t_order":2,"limit":limit,"mpm_mapping":mpm_mapping,"compute_CFL":False}
##DGMPM: Discontinous Galerkin Material Point Method
DGMPM3 = dict(parameters)
execfile('dgmpm/hyperelasticity.py', DGMPM3)


#############################################################################
#########################  Comparison  ######################################
#############################################################################

####Animated plot ###########################################################
from matplotlib import rcParams
rcParams['axes.labelsize'] = 16
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] = 16
rcParams['legend.fontsize'] = 16


if sigd>0.:
    frames=[]
    frmpm=[]
    fr2ppc=[]
    for i in MPM["time"]:
        for j in DGMPM["time"]:
            if abs(i-j)<5.e-7:
                ndg = np.where(DGMPM["time"]==j)[0][0]
                nmpm = np.where(MPM["time"]==i)[0][0]
                
                tdg = '%.2e' % j ; tmpm = '%.2e' % i
                print "dgmpm increment ",ndg
                print "mpm increment ",nmpm, " ; time difference: ",'%.2e' %abs(j-i)
                print "mpm time: ",tmpm," ; dgmpm time:",tdg
                frames.append(ndg)
                frmpm.append(nmpm)
else :
    frames=[40,80]
    frmpm=[80,160]
    start=0

#pdb.set_trace()
# frames=[5,20]
# frames=[]
# frames=[20,30,45]
subtitles=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)']

if sigd>0.:
    if sigd==5.*Sigy:
        frames=[43,89]
        frmpm=[85,176]
        #frmpm=[2*85,2*176]
        load=5
    elif sigd==50.*Sigy:
        frames=[44,88]
        frmpm=[80,160]
        #frmpm=[160,320]
        load=50
    start=0
titles=[]
for i,n1 in enumerate(frames[start:]):
    time = '%.2e' % DGMPM["time"][n1]
    mpm=frmpm[start+i]
    print n1,mpm," time t=",time

    plt.plot(MPM["pos"][:,mpm],MPM["Pi"][:,mpm],'b-x',lw=2.,ms=8.,label='MPM')
    plt.plot(DGMPM["pos"][:,n1],DGMPM["Pi"][:,n1],'g-o',lw=2.,ms=8.,label='DGMPM')
    plt.plot(DGMPM3["pos"][:,n1],DGMPM3["Pi"][:,n1],'r-',lw=2.,ms=8.,label='DGMPM RK2')
    plt.plot(DGMPM["pos"][:,n1],DGMPM["Pi_th"][:,n1],'k',lw=2.,ms=8.,label='exact')
    #plt.plot(MPM2["pos"][:,frmpm[i]],MPM2["Pi"][:,frmpm[i]],'y-x',lw=2.,ms=8.,label='MPM 2ppc')
    plt.plot(DGMPM2["pos"][:,2*n1],DGMPM2["Pi"][:,2*n1],'r-o',lw=2.,ms=8.,label='DGMPM 2ppc')
    plt.title('Contrainte longitudinale dans la barre au temps t='+str(time)+' s.',size=24.)
    plt.xlabel('x (m)',size=24.)
    plt.ylabel(r'$\sigma (Pa)$',size=28.)
    plt.legend(numpoints=1)
    plt.grid()
    #plt.show()
    legend=['mpm','dgmpm','dgmpm 2ppc','dgmpm 2ppc (RK2)','exact']
    temps=time[:-4]
    if sigd>0.:
        case=str(load)+'shock'
    else:
        case='rare'
    subtitle=subtitles[i]+r' $t = '+str(temps)+r'\times 10^{-'+str(time[-1])+'} $ s.'
    titles.append(subtitle)
    #export2DTeXFile(str(path)+'/he_stress_'+case+str(n1)+'.tex',np.array([MPM["pos"][:,mpm],DGMPM["pos"][:,n1],DGMPM2["pos"][:,2*n1],DGMPM3["pos"][:,n1],DGMPM["pos"][:,n1]]),'$x (m)$',r'$\Pi (Pa)$',str(subtitle),np.array([MPM["Pi"][:,mpm],DGMPM["Pi"][:,n1],DGMPM2["Pi"][:,2*n1],DGMPM3["Pi"][:,n1],DGMPM["Pi_th"][:,n1]]),legend)

    # export2DTeXFile(str(path)+'/dgmpm_mpm_velo'+case+str(n1)+'.tex',np.array([MPM["pos"][:,2*n1],MPM2["pos"][:,2*n1],DGMPM["pos"][:,n1],DGMPM2["pos"][:,2*n1],DGMPM3["pos"][:,n1]]),'$x (m)$','$v (m/s)$',str(subtitle),np.array([MPM["velo"][:,2*n1],MPM2["velo"][:,2*n1],DGMPM["velo"][:,n1],DGMPM2["velo"][:,2*n1],DGMPM3["velo"][:,n1]]),legend)

fileName=str(path)+'/he_stress_'+case+'.tex'
Exact=dict();Exact["pos"]=DGMPM["pos"];Exact["Pi"]=DGMPM["Pi_th"]
containers=np.array([MPM,DGMPM,DGMPM2,DGMPM3,Exact])
rowFields=['Pi']
colFields=np.array([[frmpm[0],frames[0],2*frames[0],frames[0],frames[0]],[frmpm[1],frames[1],2*frames[1]+1,frames[1],frames[1]]])
legend=['mpm','dgmpm','dgmpm 2ppc','dgmpm 2ppc (RK2)','exact']
Ylabels=[r'$\Pi (Pa)$']

#export2DGroupplot(fileName,containers,rowFields,colFields,titles,Ylabels,legend)
print DGMPM2["time"][[2*frames[0],2*frames[1]+1]],DGMPM["time"][frames]

####################################################################
fig, (ax1) = plt.subplots(1,1)

# intialize two line objects (one in each axes)
line1, = ax1.plot([], [],'r-o', ms=2.5,label='dgmpm')
line2, = ax1.plot([], [],'g--' ,ms=1.5,label='mpm')
line3, = ax1.plot([], [],'rs' , ms=2,label='dgmpm')
line4, = ax1.plot([], [],'b--s' , ms=2,label='mpm')
line5, = ax1.plot([], [],'k' , ms=2,label='fvm')
line = [line1, line2,line3,line4,line5]

ax1.grid()
ax1.set_xlabel('x (m)', fontsize=18)
ax1.set_ylabel(r'$\sigma$', fontsize=18)

ax1.legend(numpoints=1)

ax1.set_xlim(0.,length)
ax1.set_ylim(1.1*np.min(MPM["Pi"]),1.1*np.max(MPM["Pi"]))
def init():
    line[0].set_data([], [])
    line[1].set_data([], [])
    line[2].set_data([], [])
    line[3].set_data([], [])
    line[4].set_data([], [])
    return line

def animate(i):
    line[0].set_data(MPM["pos"][:,i],MPM["Pi"][:,i])
    return line

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=MPM["increments"], interval=100, blit=True)

plt.show()

