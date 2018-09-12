#!/usr/bin/pyton

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import pdb
import os
from ExactplaneWaveKin import *

directory=os.path.basename(__file__)[:-3]

if not os.path.exists('texFiles/'+str(directory)):
    os.system('mkdir texFiles/'+str(directory))
path='texFiles/'+str(directory)
"""
Comparison of the implementation of the 1D elastic-plastic set of equations
in dynamics for linear isotropic/kinematic hardening materials:
- with the MPM
- with the DGMPM
- with the FVM
"""
def export2DTeXFile(fileName,xFields,xlabel,ylabel,subtitle,yfields,*kwargs):
    TeXFile=open(fileName,"w")
    n_fields = np.shape(yfields)[0]
    n_labels = np.shape(kwargs)[0]
    # Define Paul Tol's colors (purple to red)
    ## 5 fields
    marker=['none','none','|','none','pentagone*','none','triangle*','none']
    style=['dashed','solid','solid','solid','solid','solid','dotted']
    thickness=['very thick','very thick','thick','thin','very thick','thin','thick','very thick']
    couleur=['Red','Blue','Purple','black','Yellow','black','Green','Orange']
    TeXFile.write(r'\begin{tikzpicture}[scale=0.8]');TeXFile.write('\n')
    TeXFile.write(r'\begin{axis}[xlabel='+str(xlabel)+',ylabel='+str(ylabel)+',ymajorgrids=true,xmajorgrids=true,legend pos=outer north east,title={'+subtitle+'},xmin=0.,xmax=6.]');TeXFile.write('\n')
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
    if subtitle[:3]=='(c)':
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
    marker=['none','none','+','none','pentagone*','none','triangle*','none']
    style=['dashed','solid','solid','solid','solid','solid','dotted']
    thickness=['very thick','very thick','thick','thin','very thick','thin','thick','very thick']
    couleur=['Red','Blue','Purple','black','Yellow','black','Green','Orange']
    maximum=np.zeros(row)
    minimum=np.zeros(row)
    # sum over rows (fields sigma of epsp)
    for i,field in enumerate(rowFields):
        maxim=[]
        minim=[]
        #if field=='epsp':pdb.set_trace()
        # sum over columns (t1,t2,etc.)
        for j in range(col):
            # sum over fields in plots (USL,DGMPM,etc.)
            for k in range(fields_in_plots):
                maxim.append(max(containers[k][field][:,colFields[j][k]]))
                minim.append(min(containers[k][field][:,colFields[j][k]]))
        maximum[i]=1.1*max(maxim)
        minimum[i]=1.1*min(minim)
    TeXFile=open(fileName,"w")
    # Define Paul Tol's colors (purple to red)
    TeXFile.write(r'\begin{tikzpicture}[scale=.9]');TeXFile.write('\n')
    TeXFile.write(r'\begin{groupplot}[group style={group size='+str(col)+' by '+str(row)+',');TeXFile.write('\n')
    TeXFile.write('ylabels at=edge left, yticklabels at=edge left,horizontal sep=2.ex,');TeXFile.write('\n')
    TeXFile.write('vertical sep=4ex,xticklabels at=edge bottom,xlabels at=edge bottom},');TeXFile.write('\n')
    if row==1:
        TeXFile.write(r'ymajorgrids=true,xmajorgrids=true,enlargelimits=0,xmin=0.,xmax=6.,'+str(Ylabels)+',xlabel=x (m),');TeXFile.write('\n')
    else:
        TeXFile.write(r'ymajorgrids=true,xmajorgrids=true,enlargelimits=0,xmin=0.,xmax=6.,xlabel=$x (m)$,');TeXFile.write('\n')
    TeXFile.write('axis on top,scale only axis,width=0.32\linewidth');TeXFile.write('\n')
    TeXFile.write(']');TeXFile.write('\n')
    for i,field in enumerate(rowFields): ## sum over rows
        for j in range(col):
            TeXFile.write(r'\nextgroupplot[')
            if i==0: TeXFile.write(r'title={'+str(titles[j])+'},ymin='+str(minimum[i])+',ymax='+str(maximum[i])+',')
            elif j==0: TeXFile.write(r'ylabel='+str(Ylabels[i])+',ymin='+str(minimum[i])+',ymax='+str(maximum[i])+',')
            elif j==col-1 and i==row-1: TeXFile.write(r'legend style={at={($(0.5,-0.45)+(0.cm,1cm)$)},legend columns=2},ymin='+str(minimum[i])+',ymax='+str(maximum[i]))
            else: TeXFile.write('ymin='+str(minimum[i])+',ymax='+str(maximum[i])+',')
            TeXFile.write(']');TeXFile.write('\n')
            for k in range(fields_in_plots):
                TeXFile.write(r'\addplot['+str(couleur[k])+','+str(style[k])+',mark='+str(marker[k])+','+thickness[k]+',mark size=3pt,mark repeat=2] coordinates{')
                #pdb.set_trace()
                #print field
                FIELD=containers[k][field][:,colFields[j][k]]
                xFields=containers[k]["pos"][:,colFields[j][k]]
                for l in range(len(FIELD)):
                    TeXFile.write('('+str(xFields[l])+','+str(FIELD[l])+') ')
                TeXFile.write('};\n')
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
Nelem = 50
E = 2.0e11
nu=0.3
lamb=E*nu/((1.+nu)*(1.-2.*nu))
mu=0.5*E/(1.+nu)
Sigy = 400.0e6
H = 10e9
rho = 7800.0
HEL = ((lamb+2.0*mu)/(2.0*mu))*Sigy
c=np.sqrt((lamb+2.*mu)/rho)
sigd =0.
v0=2.*HEL/(rho*c)
factor=1.
timeOut = 1.*length/c
t_order=1
timeUnload = 2*timeOut
algo = 'USL'
update_position=False
mpm_mapping=True
limit=-1
hardening='kinematic'
fvmlimiter=-1

parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"nu":nu,"Sigy":Sigy, "H":H,"rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"algo":algo,"t_order":t_order,"limit":limit,"mpm_mapping":mpm_mapping,"compute_CFL":False,"hardening":hardening,"fvmlimiter":fvmlimiter}
#################

##MPM: Material Point Method
USL = dict(parameters)
print 'Computing  USF'
execfile('mpm/planeWave.py', USL)

parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"nu":nu,"Sigy":Sigy, "H":H,"rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"algo":'USF',"t_order":t_order,"limit":limit,"mpm_mapping":mpm_mapping,"compute_CFL":False,"hardening":hardening,"fvmlimiter":fvmlimiter}
#################

##MPM: Material Point Method
USF = dict(parameters)
print 'Computing  USF'
execfile('mpm/planeWave.py', USF)


##DGMPM: Discontinuous Galerkin Material Point Method
DGMPM = dict(parameters)
print 'Computing  DGMPM (ep solver)'
execfile('dgmpm/planeWaveEP_EPsolver.py', DGMPM)

##DGMPM: Discontinuous Galerkin Material Point Method
DGMPM2 = dict(parameters)
print 'Computing  DGMPM (ac solver)'
execfile('dgmpm/planeWaveEP_ACsolver.py', DGMPM2)



# ppc=2
# #Nelem = 50/ppc
# parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy, "H":H,"rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"algo":algo,"t_order":t_order,"limit":limit,"mpm_mapping":mpm_mapping,"compute_CFL":False}

# print "=============== 2PPC COMPUTATIONS ===================="

# ##MPM: Material Point Method
# MPM2 = dict(parameters)
# print 'Computing MPM'
# execfile('mpm/elasticity.py', MPM2)

# ##DGMPM: Discontinous Galerkin Material Point Method
# DGMPM2 = dict(parameters)
# print 'Computing DGMPM'
# execfile('dgmpm/elasticity.py', DGMPM2)

# parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy, "H":H,"rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"algo":algo,"t_order":2,"limit":limit,"mpm_mapping":mpm_mapping,"compute_CFL":False}
# ##DGMPM: Discontinous Galerkin Material Point Method
# DGMPM3 = dict(parameters)
# print 'Computing DGMPM (RK2)'
# execfile('dgmpm/elasticity.py', DGMPM3)

#############################################################################
#########################  Comparison  ######################################
#############################################################################

####Animated plot ###########################################################
from matplotlib import rcParams
rcParams['axes.labelsize'] = 16
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] = 16
rcParams['legend.fontsize'] = 16


subtitles=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)']
frames=[45,65,95]
frames=[20,30,45]
titles=[]
sig_th=np.zeros((len(DGMPM["pos"][:,0]),len(frames)))
epsp_th=np.zeros((len(DGMPM["pos"][:,0]),len(frames)))

for i,n1 in enumerate(frames):
    time = '%.2e' % DGMPM["time"][n1]
    fig, (ax1, ax2) = plt.subplots(2,1)
    if DGMPM["time"][n1]<=0.5*length/c :
        temps=DGMPM["time"][n1]
    else:
        temps=DGMPM["time"][n1-1]
    
    if hardening=='isotropic':
        Sexact,Epexact,Vexact = computeAnalyticalSolutionISO(DGMPM["pos"][:,n1],length,c,temps,abs(v0),HEL,lamb,mu,H,rho)
    elif hardening=='kinematic':
        Sexact,Epexact,Vexact = computeAnalyticalSolutionKIN(DGMPM["pos"][:,n1],length,c,temps,abs(v0),HEL,lamb,mu,H,rho)
    sig_th[:,i]=-np.sign(v0)*Sexact
    epsp_th[:,i]=-np.sign(v0)*Epexact

    ax1.plot(USL["pos"][:,2*n1],USL["sig"][:,2*n1],'y',lw=2.,ms=4.,label='USL')
    ax1.plot(USF["pos"][:,2*n1],USF["sig"][:,2*n1],'b',lw=2.,ms=4.,label='USF')
    ax1.plot(DGMPM["pos"][:,n1],DGMPM["sig"][:,n1],'r',lw=2.,ms=4.,label='DGMPM (ep solver)')
    ax1.plot(DGMPM2["pos"][:,n1],DGMPM2["sig"][:,n1],'r--',lw=2.,ms=4.,label='DGMPM (ac solver)')
    ax1.plot(DGMPM["pos"][:,n1],-np.sign(v0)*Sexact,'k',lw=2.,ms=4.,label='exact')
    ax1.set_title('Stress at time t='+str(time)+' s.',size=24.)
    ax1.set_xlabel('x (m)')
    ax1.set_ylabel(r'$\sigma (Pa)$')

    ax2.plot(USL["pos"][:,2*n1],USL["epsp"][:,2*n1],'y',lw=2.,ms=4.,label='USL')
    ax2.plot(USF["pos"][:,2*n1],USF["epsp"][:,2*n1],'b',lw=2.,ms=4.,label='USF')
    ax2.plot(DGMPM["pos"][:,n1],DGMPM["epsp"][:,n1],'r',lw=2.,ms=4.,label='DGMPM (ep solver)')
    ax2.plot(DGMPM2["pos"][:,n1],DGMPM2["epsp"][:,n1],'r--',lw=2.,ms=4.,label='DGMPM (ac solver)')
    ax2.plot(DGMPM["pos"][:,n1],-np.sign(v0)*Epexact,'k',lw=2.,ms=4.,label='exact')
    ax2.set_title('Plastic strain at time t='+str(time)+' s.')
    ax2.set_xlabel('x (m)')
    ax2.set_ylabel(r'$\varepsilon^p (Pa)$')
    ax1.legend(numpoints=1)
    ax1.grid();ax2.grid()
    #plt.show()
    legend=['usl 1ppc','dgmpm 1ppc (ep solver)','dgmpm 1ppc (ac solver)','exact']
    temps=time[:-4]
    subtitle=subtitles[i]+r' $t = '+str(temps)+r'\times 10^{-'+str(time[-1])+'} $ s.'
    titles.append(subtitle)
    #export2DTeXFile(str(path)+'/EP_dgmpm_mpm_stress'+str(n1)+'.tex',np.array([USL["pos"][:,2*n1],DGMPM["pos"][:,n1],DGMPM2["pos"][:,n1],DGMPM["pos"][:,n1]]),'$x (m)$',r'$\sigma (Pa)$',str(subtitle),np.array([USL["sig"][:,2*n1],DGMPM["sig"][:,n1],DGMPM2["sig"][:,n1],-np.sign(v0)*Sexact]),legend)
    #export2DTeXFile(str(path)+'/EP_dgmpm_mpm_epsp'+str(n1)+'.tex',np.array([USL["pos"][:,2*n1],DGMPM["pos"][:,n1],DGMPM2["pos"][:,n1],DGMPM["pos"][:,n1]]),'$x (m)$',r'$\eps^p$',str(subtitle),np.array([USL["epsp"][:,2*n1],DGMPM["epsp"][:,n1],DGMPM2["epsp"][:,n1],-np.sign(v0)*Epexact]),legend)
    #export2DTeXFile(str(path)+'/EP_dgmpm_mpm_velo'+str(n1)+'.tex',np.array([USL["pos"][:,2*n1],DGMPM["pos"][:,n1],DGMPM2["pos"][:,n1],DGMPM["pos"][:,n1]]),'$x (m)$',r'$\eps^p$',str(subtitle),np.array([USL["velo"][:,2*n1],DGMPM["Velocity"][:,n1],DGMPM2["Velocity"][:,n1],-np.sign(v0)*Vexact]),legend)

fileName=str(path)+'/ep_dgmpm_mpm.tex'
Exact=dict();Exact["pos"]=DGMPM["pos"];Exact["sig"]=sig_th;Exact["epsp"]=epsp_th
# USL["pos"][:,2*n1],DGMPM["pos"][:,n1],DGMPM2["pos"][:,n1],DGMPM["pos"][:,n1]
containers=np.array([USL,DGMPM,DGMPM2,Exact])
rowFields=['sig','epsp']
colFields=np.array([[40,20,20,0],[60,30,30,1],[90,45,45,2]])
legend=['usl 1ppc','dgmpm 1ppc (elastic-plastic solver)','dgmpm 1ppc (elastic solver)','exact']
Ylabels=[r'$\sigma (Pa)$',r'$\eps^p $']

export2DGroupplot(fileName,containers,rowFields,colFields,titles,Ylabels,legend)

    
"""
####################################################################
fig, (ax1, ax2) = plt.subplots(2,1)

# intialize two line objects (one in each axes)
line1, = ax1.plot([], [],'b-o' ,ms=1.5,label='ep')
line2, = ax1.plot([], [],'k' , ms=2,label='fvm')
line3, = ax1.plot([], [],'k-o' , ms=2,label='fvm SB')
line4, = ax1.plot([], [],'g-+' , ms=2,label='fem')

line5, = ax2.plot([], [],'b-o' , ms=2,label='ep')
line6, = ax2.plot([], [],'k', ms=2.5,label='fvm')
line7, = ax2.plot([], [],'k-o', ms=2.5,label='fvm SB')
line8, = ax2.plot([], [],'g-+' , ms=2,label='fem')


line = [line1, line2,line3,line4,line5,line6,line7,line8]

ax1.grid()
ax2.grid()
ax1.set_xlabel('x (m)', fontsize=18)
ax2.set_xlabel('x (m)', fontsize=18)
ax2.set_ylabel(r'$\varepsilon^p$', fontsize=18)
ax1.set_ylabel(r'$\sigma$', fontsize=18)

ax1.set_xlim(0.,length)
ax2.set_xlim(0.,length)
ax1.set_ylim(1.1*np.min(FEM["sig"]),1.1*np.max(FEM["sig"]))
ax2.set_ylim(1.1*np.min(FEM["epsp"]),1.1*np.max(FEM["epsp"]))
ax1.legend(numpoints=1)
def init():
    line[0].set_data([], [])
    line[1].set_data([], [])
    line[2].set_data([], [])
    line[3].set_data([], [])
    line[4].set_data([], [])
    line[5].set_data([], [])
    line[6].set_data([], [])
    line[7].set_data([], [])
    return line

def animate(i):
    line[0].set_data(DGMPM["pos"][:,i],DGMPM["sig"][:,i])
    line[1].set_data(FVM["centroids"],FVM["sig"][:,i])
    line[2].set_data(FVM2["centroids"],FVM2["sig"][:,i])
    line[3].set_data(FEM["centroids"],FEM["sig"][:,i])

    line[4].set_data(DGMPM["pos"][:,i],DGMPM["epsp"][:,i])
    line[5].set_data(FVM["centroids"],FVM["EP"][:,i])
    line[6].set_data(FVM2["centroids"],FVM2["EP"][:,i])
    line[7].set_data(FEM["centroids"],FEM["epsp"][:,i])
    return line

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=DGMPM["increments"], interval=100, blit=True)

plt.show()
"""
