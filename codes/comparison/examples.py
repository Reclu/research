
#!/usr/bin/pyton

import numpy as np
import matplotlib.pyplot as plt
from exportToTeXFile import *

x=np.array([np.linspace(2,10,50),np.linspace(10,20,50)])
x=x.T

sqrt=np.sqrt(x)
log=np.log(x)
cos=5.*np.cos(x)
sin=5.*np.sin(x)

## Exporting one TeX file
legend=[r'$f=\sqrt{x}$',r'$f=\log{x}$',r'$f=\cos{x}$',r'$f=\sin{x}$']
filename='example_export2DTeXFile.tex'
export2DTeXFile(filename,np.array([x[:,0],x[:,0],x[:,0],x[:,0]]),'$x (m)$',r'Functions','',np.array([sqrt[:,0],log[:,0],cos[:,0],sin[:,0]]),legend)


## Exporting groupplots: first line for f(x) and second line for f(2x); first column for 2<x<10, and second column for 10<x<20
sqrt2=np.sqrt(2.*x)
log2=np.log(2.*x)
cos2=5.*np.cos(2.*x)
sin2=5.*np.sin(2.*x)


SQRT=dict();
SQRT["x"]=x
SQRT["f1"]=sqrt
SQRT["f2"]=sqrt2

LOG=dict();
LOG["x"]=x
LOG["f1"]=log
LOG["f2"]=log2

COS=dict();
COS["x"]=x
COS["f1"]=cos
COS["f2"]=cos2


SIN=dict();
SIN["x"]=x
SIN["f1"]=sin
SIN["f2"]=sin2

containers=np.array([SQRT,LOG,COS,SIN])
filename='example_export2DGroupplot.tex'
rowFields=['f1','f2']
colFields=np.array([[0,0,0,0],[1,1,1,1]])
titles=['(a)$ 2<x<10$ ','(b) $10<x<20$']
legend=[r'$f=\sqrt{}$',r'$f=\log{}$',r'$f=\cos{}$',r'$f=\sin{}$']
Ylabels=[r'$f(x)$',r'$f(2x)$']
Xlabels='$x$'

export2DGroupplot(filename,containers,rowFields,colFields,titles,Ylabels,Xlabels,legend)
