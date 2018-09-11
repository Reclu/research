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
- with the MPM
- with the DGMPM
"""
def export2DTeXFile(fileName,xFields,xlabel,ylabel,subtitle,yfields,*kwargs):
    TeXFile=open(fileName,"w")
    n_fields = np.shape(yfields)[0]
    n_labels = np.shape(kwargs)[0]
    # Define Paul Tol's colors (purple to red)
    marker=['none','none','*','none','|','x','pentagone*','none','triangle*']
    style=['dashed','densely dotted','solid','solid','solid','only marks','solid','solid']
    thickness=['very thick','very thick','thick','very thick','very thick','thick','thin','thin','thick']
    couleur=['Red','Orange','Duck','Blue','Purple','Green','black','Yellow','black','Green']
    TeXFile.write(r'\begin{tikzpicture}[scale=1.]');TeXFile.write('\n')
    TeXFile.write(r'\begin{axis}[xlabel='+str(xlabel)+',ylabel='+str(ylabel)+',ymajorgrids=true,xmajorgrids=true,legend pos=outer north east,title={'+subtitle+'}]');TeXFile.write('\n')
    legend=''
    for i in range(n_fields):
        if i==0:
            legend=legend+kwargs[0][i]
        else:
            legend=legend+','+kwargs[0][i]
        TeXFile.write(r'\addplot['+str(couleur[i])+','+str(thickness[i])+',mark='+str(marker[i])+','+str(style[i])+',mark size=3pt] coordinates {')
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
    TeXFile=open(fileName,"w")
#     n_fields = np.shape(yFields)[0];
#     n_labels = np.shape(kwargs)[0]
    # Define Paul Tol's colors (purple to red)
    marker=['none','none','*','none','|','x','pentagone*','none','triangle*']
    style=['dashed','densely dotted','solid','solid','solid','only marks','solid','solid']
    thickness=['very thick','very thick','thick','very thick','very thick','thick','thin','thin','thick']
    couleur=['Red','Orange','Duck','Blue','Purple','Green','black','Yellow','black','Green']
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
    TeXFile.write(r'\begin{tikzpicture}[scale=.9]');TeXFile.write('\n')
    TeXFile.write(r'\begin{groupplot}[group style={group size='+str(col)+' by '+str(row)+',');TeXFile.write('\n')
    TeXFile.write('ylabels at=edge left, yticklabels at=edge left,horizontal sep=4.ex,');TeXFile.write('\n')
    TeXFile.write('vertical sep=2ex,xticklabels at=edge bottom,xlabels at=edge bottom},');TeXFile.write('\n')
    if row==1:
        TeXFile.write(r'ymajorgrids=true,xmajorgrids=true,enlargelimits=0,xmin=0.,xmax=6.,'+str(Ylabels)+',xlabel=x (m),');TeXFile.write('\n')
    else:
        TeXFile.write(r'ymajorgrids=true,xmajorgrids=true,enlargelimits=0,xmin=0.,xmax=6.,xlabel=x (m),');TeXFile.write('\n')
    TeXFile.write('axis on top,scale only axis,width=0.45\linewidth');TeXFile.write('\n')
    TeXFile.write(']');TeXFile.write('\n')
    for i,field in enumerate(rowFields): ## sum over rows
        for j in range(col):
            TeXFile.write(r'\nextgroupplot[')
            if i==0: TeXFile.write(r'title={'+str(titles[j])+'},ymin='+str(minimum[i])+',ymax='+str(maximum[i])+',')
            elif j==0: TeXFile.write(r'ylabel='+str(Ylabels[i])+',ymin='+str(minimum[i])+',ymax='+str(maximum[i])+',')
            elif j==col-1 and i==row-1: TeXFile.write(r'legend style={at={($(0.62,-0.35)+(0.9cm,1cm)$)},legend columns=4},ymin='+str(minimum[i])+',ymax='+str(maximum[i]))
            else: TeXFile.write('ymin='+str(minimum[i])+',ymax='+str(maximum[i])+',')
            TeXFile.write(']');TeXFile.write('\n')
            for k in range(fields_in_plots):
                TeXFile.write(r'\addplot['+str(couleur[k])+','+str(style[k])+',mark='+str(marker[k])+','+thickness[k]+',mark size=3pt,mark repeat=2] coordinates{')
                #pdb.set_trace()
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
Sigy = 400.0e6
H = 10e9
rho = 7800.0
c=np.sqrt(E/rho)
sigd =0.
v0=0.5*Sigy/(2*rho*c)
factor=1.
timeOut = 0.5*length/(np.sqrt(E/rho))
t_order=1
timeUnload = 2*timeOut
algo = 'USL'
update_position=False
mpm_mapping=True
limit=-1
hardening='isotropic'


parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy, "H":H,"rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"algo":algo,"t_order":t_order,"limit":limit,"mpm_mapping":mpm_mapping,"compute_CFL":False,"hardening":hardening}
#################


##MPM: Material Point Method
MPM = dict(parameters)
print 'Computing MPM'
execfile('mpm/elasticity.py', MPM)

mpm_mapping=False
parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy, "H":H,"rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"algo":algo,"t_order":t_order,"limit":limit,"mpm_mapping":mpm_mapping,"compute_CFL":False,"hardening":hardening}
#################


##MPM: Material Point Method
PIC = dict(parameters)
print 'Computing MPM'
execfile('mpm/elasticity.py', PIC)

mpm_mapping=True
parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy, "H":H,"rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"algo":algo,"t_order":t_order,"limit":limit,"mpm_mapping":mpm_mapping,"compute_CFL":False,"hardening":hardening}

##DGMPM: Discontinuous Galerkin Material Point Method
DGMPM = dict(parameters)
print 'Computing  DGMPM'
execfile('dgmpm/elasticity.py', DGMPM)


ppc=2
#Nelem = 50/ppc
parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy, "H":H,"rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"algo":algo,"t_order":t_order,"limit":limit,"mpm_mapping":mpm_mapping,"compute_CFL":False,"hardening":hardening}

print "=============== 2PPC COMPUTATIONS ===================="

##MPM: Material Point Method
MPM2 = dict(parameters)
print 'Computing MPM'
execfile('mpm/elasticity.py', MPM2)

#Nelem = 50/ppc
parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy, "H":H,"rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"algo":algo,"t_order":t_order,"limit":limit,"mpm_mapping":False,"compute_CFL":False,"hardening":hardening}

print "=============== 2PPC COMPUTATIONS ===================="

##MPM: Material Point Method
PIC = dict(parameters)
print 'Computing MPM (PIC)'
execfile('mpm/elasticity.py', PIC)


##DGMPM: Discontinous Galerkin Material Point Method
DGMPM2 = dict(parameters)
print 'Computing DGMPM'
execfile('dgmpm/elasticity.py', DGMPM2)

parameters = {"CFL":CFL,"Nelem":Nelem,"NTmaxi":NTmaxi,"ppc":ppc,"length":length,"Young":E,"Sigy":Sigy, "H":H,"rho":rho,"sigd":sigd,"timeOut":timeOut,"timeUnload":timeUnload,"update_position":update_position,"v0":v0,"factor":factor,"algo":algo,"t_order":2,"limit":limit,"mpm_mapping":mpm_mapping,"compute_CFL":False,"hardening":hardening}
##DGMPM: Discontinous Galerkin Material Point Method
DGMPM3 = dict(parameters)
print 'Computing DGMPM (RK2)'
execfile('dgmpm/elasticity.py', DGMPM3)



#############################################################################
#########################  Comparison  ######################################
#############################################################################

####Animated plot ###########################################################
from matplotlib import rcParams
rcParams['axes.labelsize'] = 16
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] = 16
rcParams['legend.fontsize'] = 16


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
print factor_dgmpms,factor_mpms
factor2=DGMPM2["CFL"]/CFL


### Energy plots
plt.plot(MPM["time"][:-1],MPM["NRG"][:-1]/max(MPM["NRG"][:-1]),'b-x',lw=2.,label='mpm 1ppc')
plt.plot(PIC["time"][:-1],PIC["NRG"][:-1]/max(PIC["NRG"][:-1]),'b-x',lw=2.,label='pic 1ppc')
plt.plot(MPM2["time"][:-1],MPM2["NRG"][:-1]/max(MPM2["NRG"][:-1]),'r-x',lw=2.,label='mpm 2ppc')
plt.plot(DGMPM["time"][:-1],DGMPM["NRG"][:-1]/max(DGMPM["NRG"][:-1]),'bo',lw=2.,label='dgmpm 1ppc')
plt.plot(DGMPM2["time"][:-1],DGMPM2["NRG"][:-1]/max(DGMPM2["NRG"][:-1]),'ro',lw=2.,label='dgmpm 2ppc')
plt.grid()
plt.legend(numpoints=1)
plt.show()

# export2pgfPlot('NRG_mpm_1ppc.pgf',MPM["time"][:-1],MPM["NRG"][:-1],'t','NRG')
# export2pgfPlot('NRG_mpm_2ppc.pgf',MPM2["time"][:-1],MPM2["NRG"][:-1],'t','NRG')
# export2pgfPlot('NRG_modmpm_1ppc.pgf',DGMPM["time"][:-1],DGMPM["NRG"][:-1],'t','NRG')
# export2pgfPlot('NRG_modmpm_2ppc.pgf',DGMPM2["time"][:-1],DGMPM2["NRG"][:-1],'t','NRG')
legend=['usl 1ppc','usl 2ppc','usl-pic 2ppc','dgmpm 1ppc','dgmpm 2ppc','dgmpm 2ppc (RK2)','exact']
export2DTeXFile(str(path)+'/dgmpm_mpm_energies.tex',np.array([MPM["time"][:-1],MPM2["time"][:-1],PIC["time"][:-1],DGMPM["time"][:-1],DGMPM2["time"][:-1],DGMPM3["time"][:-1],np.array([0.,1.e-8])]),'$time (s)$',r'$\frac{e}{e_{max}}$','',np.array([MPM["NRG"][:-1]/max(MPM["NRG"][:-1]),MPM2["NRG"][:-1]/max(MPM2["NRG"][:-1]),PIC["NRG"][:-1]/max(PIC["NRG"][:-1]),DGMPM["NRG"][:-1]/max(DGMPM["NRG"][:-1]),DGMPM2["NRG"][:-1]/max(DGMPM2["NRG"][:-1]),DGMPM3["NRG"][:-1]/max(DGMPM3["NRG"][:-1]),np.array([1.,1.])]),legend)
titles=[]
frames=[5,20]
for n1 in frames:
    time = '%.2e' % MPM["time"][2*n1]
    plt.plot(MPM["pos"][:,n1],MPM["Sth"][:,2*n1],'k-',lw=2.,ms=8.,label='analytical')
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
    legend=['usl 1ppc','usl-pic 1ppc','usl 2ppc','dgmpm 1ppc','dgmpm 2ppc','dgmpm 2ppc (RK2)','exact']
    temps=time[:-4]
    if n1==5 :
        subtitle=r'(a) time $t = '+str(temps)+r'\times 10^{-'+str(time[-1])+'} $ s.'
        titles.append(subtitle)

    if n1==20 :
        subtitle=r'(b) time $t = '+str(temps)+r'\times 10^{-'+str(time[-1])+'} $ s.'
        titles.append(subtitle)

    #export2DTeXFile(str(path)+'/dgmpm_mpm_diffusion'+str(n1)+'.tex',np.array([MPM["pos"][:,2*n1],MPM2["pos"][:,2*n1],DGMPM["pos"][:,n1],DGMPM2["pos"][:,2*n1],DGMPM3["pos"][:,n1],MPM["pos"][:,2*n1]]),'$x (m)$',r'$\sigma (Pa)$',str(subtitle),np.array([MPM["sig"][:,2*n1],MPM2["sig"][:,2*n1],DGMPM["sig"][:,n1],DGMPM2["sig"][:,2*n1],DGMPM3["sig"][:,n1],MPM["Sth"][:,2*n1]]),legend)
    #export2DTeXFile(str(path)+'/dgmpm_mpm_velo'+str(n1)+'.tex',np.array([MPM["pos"][:,2*n1],MPM2["pos"][:,2*n1],DGMPM["pos"][:,n1],DGMPM2["pos"][:,2*n1],DGMPM3["pos"][:,n1],DGMPM2["pos"][:,2*n1]]),'$x (m)$','$v (m/s)$',str(subtitle),np.array([MPM["velo"][:,2*n1],MPM2["velo"][:,2*n1],DGMPM["velo"][:,n1],DGMPM2["velo"][:,2*n1],DGMPM3["velo"][:,n1],DGMPM2["Vth"][:,2*n1]]),legend)

fileName=str(path)+'/dgmpm_mpm_elasticity.tex'

Exact=dict();Exact["pos"]=DGMPM["pos"];Exact["sig"]=DGMPM["Sth"];Exact["velo"]=DGMPM["Vth"]
containers=np.array([MPM,MPM2,PIC,DGMPM,DGMPM2,DGMPM3,Exact])
rowFields=['sig','velo']
colFields=np.array([[10,10,10,5,10,5,5],[20,20,20,10,20,10,10]])
#titles=['(a)','(b)']
legend=['usl 1ppc','usl 2ppc','usl-pic 2ppc','dgmpm 1ppc','dgmpm 2ppc','dgmpm 2ppc (RK2)','exact']
Ylabels=[r'$\sigma (Pa)$','v (m/s)']

export2DGroupplot(fileName,containers,rowFields,colFields,titles,Ylabels,legend)

