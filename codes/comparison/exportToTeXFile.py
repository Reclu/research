#!/usr/bin/pyton

import numpy as np
import pdb


"""
This files gathers data export methods that can be used within a LaTeX file providing the use of:
\usepackage{pgfplots} 
pgfplot is a very powerfull LaTeX package for scientific plotting (see http://sourceforge.net/projects/pgfplots)

"""

def export2DTeXFile(fileName,xFields,xlabel,ylabel,subtitle,yFields,*kwargs):
    """
    Writing of .tex file using the environment axis embedded in a tikzpicture environment.
    Given numpy arrays containing points coordinates (X,Y), several curves can be plotted for comparison purposes.
    Hence, different colors, styles and markers are set for each curve.

    xFields is an array containing X values for each curve. 
    Those values are not necessarily the same for all the curves plotted (it is for instance the case if one compares numerical methods based on different space discretization so that the positions of nodes are not the same i.e: finite element nodes and material points)

    xFields then stands for Y values.

    The axis names are given by xlabel and ylabel parameteres and the legend is given by the list *kwargs

    A subtitle can also be used.

    The files make use of Paul Tol's colors (see https://personal.sron.nl/~pault/)
    """
    TeXFile=open(fileName,"w")
    n_fields = np.shape(yFields)[0]
    n_labels = np.shape(kwargs)[0]
    # Define Paul Tol's colors (purple to red)
    marker=['none','none','*','none','|','x','pentagone*','none','triangle*']
    style=['dashed','densely dotted','solid','solid','solid','only marks','solid','solid']
    thickness=['very thick','very thick','thick','very thick','very thick','thick','thin','thin','thick']
    couleur=['Red','Orange','Duck','Blue','Purple','Green','black','Yellow','black','Green']
    TeXFile.write(r'\begin{tikzpicture}[scale=1.]');TeXFile.write('\n')
    TeXFile.write(r'\begin{axis}[xlabel='+str(xlabel)+',ylabel='+str(ylabel)+',ymajorgrids=true,xmajorgrids=true,legend pos=outer north east,title={'+subtitle+'}]');
    TeXFile.write('\n')
    legend=''
    for i in range(n_fields):
        if i==0:
            legend=legend+kwargs[0][i]
        else:
            legend=legend+','+kwargs[0][i]
        TeXFile.write(r'\addplot['+str(couleur[i])+','+str(thickness[i])+',mark='+str(marker[i])+','+str(style[i])+',mark size=3pt] coordinates {')
        for j in range(np.shape(yFields[i])[0]):
            TeXFile.write('('+str(xFields[i][j])+','+str(yFields[i][j])+') ')
        TeXFile.write('};\n')
    if subtitle[:3]=='(c)':
        TeXFile.write(r'\legend{'+str(legend)+'}')
    else:
        TeXFile.write(r'%\legend{'+str(legend)+'}')
    TeXFile.write('\n')    
    TeXFile.write(r'\end{axis}')
    TeXFile.write('\n')
    TeXFile.write('\end{tikzpicture}')
    TeXFile.close()


def export2DGroupplot(fileName,containers,rowFields,colFields,titles,Ylabels,Xlabels,legend):
    """
    Writing of .tex file using the environment groupplot embedded in a tikzpicture environment.
    \usetikzlibrary{pgfplots.groupplots} is thus required in the main LaTeX file

    The process is similar to export2DTeXFile() but here, many graphs are plotted.

    Hence, use is made of dict() that are gathered in containers.
    It is assumed that those dictionnaries contain several fields arrays/matrices (position, velocity, stress etc.)
    e.g. : each column of the array "position" is related to one time step and each row is related to one node.
    
    Then rowFields is a list containing the keys that allows to get the y-axis field for each line of the groupplot.
    For instance, consider a groupplot for which 4 numerical methods are compared in terms of velocity and pressure along a line for two 3 time steps.
    First, containers contains 4 dicts, one for each method. 
    Second, considering a 2-rows x 3-columns groupplot, that is:
    
    ----------    ----------   ----------                                                                      
    |        |    |        |   |        |
    |  v  t1 |    |  v  t2 |   |  v  t3 |
    |        |    |        |   |        |  
    ----------    ----------   ---------- 
    ----------    ----------   ----------                                                                      
    |        |    |        |   |        |
    |  p  t1 |    |  p  t2 |   |  p  t3 |
    |        |    |        |   |        |  
    ----------    ----------   ---------- 
    
    rowFields is ['velocity','pressure']
    At last, colFields contain the indices corresponding to the times t1 t2 and t3 within the arrays, velocity and pressure, namely
    if t1 corresponds to increment 20 in all numerical scheme exept one for which it is 25; t2 is 50 or 55; and t3 is 60 and 65,
    colFields=np.array([[20,20,20,25],[50,50,50,55],[60,60,60,65]])

    The legend is also a list of label holding for each curve.
    """
    row=len(rowFields)
    col=len(colFields)
    fields_in_plots=len(containers)
    TeXFile=open(fileName,"w")
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
    TeXFile.write(r'\begin{groupplot}[group style={group size='+str(col)+' by '+str(row)+',')
    TeXFile.write('\n')
    TeXFile.write('ylabels at=edge left, yticklabels at=edge left,horizontal sep=4.ex,')
    TeXFile.write('\n')
    TeXFile.write('vertical sep=2ex,xticklabels at=edge bottom,xlabels at=edge bottom},')
    TeXFile.write('\n')
    if row==1:
        TeXFile.write(r'ymajorgrids=true,xmajorgrids=true,enlargelimits=0,xmin=0.,xmax=6.,ylabel='+str(Ylabels[0])+',xlabel='+str(Xlabel)+',')
        TeXFile.write('\n')
    else:
        TeXFile.write(r'ymajorgrids=true,xmajorgrids=true,enlargelimits=0,xmin=0.,xmax=6.,xlabel='+str(Xlabel)+',')
        TeXFile.write('\n')
    TeXFile.write('axis on top,scale only axis,width=0.45\linewidth')
    TeXFile.write('\n')
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
        TeXFile.write(r'\addlegendentry{'+str(lab)+'}')
        TeXFile.write('\n')
    TeXFile.write('\n')    
    TeXFile.write(r'\end{groupplot}')
    TeXFile.write('\n')
    TeXFile.write('\end{tikzpicture}')
    TeXFile.close()


def export2pgfPlot(fileName,xfield,yfield,xlabel,ylabel):
    """
    Writing of a table file containing point values that can be included within pgfplots through the commands "table" or "file"
    """
    dataFile=open(fileName,"w")
    dataFile.write('# Curve ('+str(xlabel)+';'+str(ylabel)+') '+str(len(xfield))+' points.\n')
    for i,x in enumerate(xfield):
        dataFile.write(str(x)+' '+str(yfield[i])+' i\n')
    dataFile.close()
