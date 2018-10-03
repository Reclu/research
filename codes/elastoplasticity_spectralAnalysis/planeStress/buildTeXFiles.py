#!/usr/bin/python
import numpy as np
import pdb
import os
import sys

def buildTeXFiles(names,pgfFiles,xlabels,ylabels,zlabels,subtitle,srcX,srcY):
    # srcX and srcY contain key to get the correct column (i.e. fields) in pgfFiles for instance "sigma_11"
    for i,nom in enumerate(names):
        TeXFile=open(nom,"w")
        ## For regular plots (i.e. not daviatoric plane)
	marker=['+','x','star','asterisk','none','none','star','pentagone*']
        style=['dashed','solid','solid','solid','solid','dashed','solid','pentagone*']
        thickness=['very thick','very thick','very thick','thick','thin','very thick','very thick','thin','thin','thick']
        couleur=['Red','Blue','Orange','Purple','Green','Duck']
        TeXFile.write(r'\begin{tikzpicture}[scale=0.9]');TeXFile.write('\n')
        if subtitle[i][:3]=='(c)':
	    TeXFile.write(r'\tikzset{cross/.style={cross out, draw=black, minimum size=2*(#1-\pgflinewidth), inner sep=0pt, outer sep=0pt},cross/.default={2.5pt}}');TeXFile.write('\n')
	    
	
            TeXFile.write(r'\begin{axis}[width=.75\textwidth,view={135}{35.2643},xlabel='+str(xlabels[i])+',ylabel='+str(ylabels[i])+',zlabel='+str(zlabels[i])+',xmin=-1.e8,xmax=1.e8,ymin=-1.e8,ymax=1.e8,axis equal,axis lines=center,axis on top,xtick=\empty,ytick=\empty,ztick=\empty,every axis y label/.style={at={(rel axis cs:0.,.5,-0.65)}, anchor=west}, every axis x label/.style={at={(rel axis cs:0.5,.,-0.65)}, anchor=east}, every axis z label/.style={at={(rel axis cs:0.,.0,.18)}, anchor=north}]');TeXFile.write('\n')
	    TeXFile.write(r'\node[below] at (1.1e8,0.,0.) {$\sigma^y$};');TeXFile.write('\n')
	    TeXFile.write(r'\node[above] at (-1.1e8,0.,0.) {$-\sigma^y$};');TeXFile.write('\n')
	    TeXFile.write(r'\draw (1.e8,0.,0.) node[cross,rotate=10] {};');TeXFile.write('\n')
	    TeXFile.write(r'\draw (-1.e8,0.,0.) node[cross,rotate=10] {};');TeXFile.write('\n')
	    TeXFile.write(r'\node[white]  at (0,0.,1.42e8) {};');TeXFile.write('\n')
	    for j,name in enumerate(pgfFiles[i][:len(couleur)]):
                #pdb.set_trace()
                if name[25:25+12]=='CylindreDevP': ##  yield surface
                    TeXFile.write(r'\addplot3+[gray,dashed,thin,no markers] file {chapter5/pgfFigures/'+name+'};')
                    TeXFile.write('\n')
                else:
                    TeXFile.write(r'\addplot3+['+couleur[j]+',very thick,no markers] file {chapter5/pgfFigures/'+name+'};')
                    #TeXFile.write(r'\addplot3+[black,thick,no markers] file {chapter5/pgfFigures/'+name+'};')
                    TeXFile.write('\n')
        else:
            if subtitle[i][:3]=='(b)':
                TeXFile.write(r'\begin{axis}[colorbar,colorbar style={title= {$p$}},xlabel='+str(xlabels[i])+',ylabel='+str(ylabels[i])+',ymajorgrids=true,xmajorgrids=true]');TeXFile.write('\n')
            else :
                TeXFile.write(r'\begin{axis}[xlabel='+str(xlabels[i])+',ylabel='+str(ylabels[i])+',ymajorgrids=true,xmajorgrids=true]');TeXFile.write('\n')
            for j,name in enumerate(pgfFiles[i]):
                #pdb.set_trace()
                if name[25:25+12]=='DPslow_yield': ##  yield surface
                    TeXFile.write(r'\addplot[gray,thin] table[x='+str(srcX[i])+',y='+str(srcY[i])+'] {chapter5/pgfFigures/'+name+'};')
                    TeXFile.write('\n')
                else:
                    TeXFile.write(r'\addplot[mesh,point meta = \thisrow{p},very thick,no markers] table[x='+str(srcX[i])+',y='+str(srcY[i])+'] {chapter5/pgfFigures/'+name+'};')
                    TeXFile.write('\n')
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

def buildTeXFiles2(names,pgfFiles,xlabels,ylabels,zlabels,srcX,srcY,ylim):
    # srcX and srcY contain key to get the correct column (i.e. fields) in pgfFiles for instance "sigma_11"
    for i,nom in enumerate(names):
        ## For regular plots (i.e. not daviatoric plane)
        marker=['star','asterisk','+','x','none','none','none','none']
        style=['dashed','solid','solid','solid','solid','dashed','solid','pentagone*']
        thickness=['very thick','very thick','very thick','thick','thin','very thick','very thick','thin','thin','thick']
        couleur=['Red','Blue','Orange','Purple','Green','Duck']
        if len(nom)!=2:
	    TeXFile=open(nom,"w")
            TeXFile.write(r'\tikzset{cross/.style={cross out, draw=black, minimum size=2*(#1-\pgflinewidth), inner sep=0pt, outer sep=0pt},cross/.default={2.5pt}}');TeXFile.write('\n')
	    TeXFile.write(r'\begin{tikzpicture}[scale=0.9]');TeXFile.write('\n')
            TeXFile.write(r'\begin{axis}[width=.75\textwidth,view={135}{35.2643},xlabel='+str(xlabels[i])+',ylabel='+str(ylabels[i])+',zlabel='+str(zlabels[i])+',xmin=-1.e8,xmax=1.e8,ymin=-1.e8,ymax=1.e8,axis equal,axis lines=center,axis on top,xtick=\empty,ytick=\empty,ztick=\empty,every axis y label/.style={at={(rel axis cs:0.,.5,-0.65)}, anchor=west}, every axis x label/.style={at={(rel axis cs:0.5,.,-0.65)}, anchor=east}, every axis z label/.style={at={(rel axis cs:0.,.0,.18)}, anchor=north},legend style={at={(.225,.59)}}]');TeXFile.write('\n')
	    TeXFile.write(r'\node[below] at (1.1e8,0.,0.) {$\sigma^y$};');TeXFile.write('\n')
	    TeXFile.write(r'\node[above] at (-1.1e8,0.,0.) {$-\sigma^y$};');TeXFile.write('\n')
	    TeXFile.write(r'\draw (1.e8,0.,0.) node[cross,rotate=10] {};');TeXFile.write('\n')
	    TeXFile.write(r'\draw (-1.e8,0.,0.) node[cross,rotate=10] {};');TeXFile.write('\n')
	    TeXFile.write(r'\node[white]  at (0,0.,1.42e8) {};');TeXFile.write('\n')
	    #TeXFile.write(r'\begin{axis}[width=.75\textwidth,view={135}{35.2643},xlabel='+str(xlabels[i])+',ylabel='+str(ylabels[i])+',zlabel='+str(zlabels[i])+',xmin=-1.e8,xmax=1.e8,ymin=-1.e8,ymax=1.e8,axis equal,axis lines=center,axis on top,ztick=\empty,legend style={at={(.225,.59)}}]');TeXFile.write('\n')
            for j,name in enumerate(pgfFiles[i][:len(couleur)]):
	    	#pdb.set_trace()
	    	if name[27:27+12]=='CylindreDevP': ##  yield surface
                    TeXFile.write(r'\addplot3+[gray,dashed,thin,no markers] file {chapter5/pgfFigures/'+name+'};')
                    TeXFile.write(r'\addlegendentry{initial yield surface}')
                    TeXFile.write('\n')
                else:
                    TeXFile.write(r'\addplot3+['+couleur[j]+',mark='+marker[j]+',mark repeat=20,mark size=3pt,very thick] file {chapter5/pgfFigures/'+name+'};\n')
		    TeXFile.write(r'\addlegendentry{loading path '+str(j+1)+'}')
                    TeXFile.write('\n')
            TeXFile.write(r'\end{axis}')
        else:
            TeXFile=open(nom[0][:11]+nom[0][-5:],"w")
            TeXFile.write(r'\begin{tikzpicture}[scale=0.9]');TeXFile.write('\n')
            TeXFile.write(r'\begin{groupplot}[group style={group size=2 by 1,');TeXFile.write('\n')
            TeXFile.write('ylabels at=edge left, yticklabels at=edge left,horizontal sep=3.ex,');TeXFile.write('\n')
            TeXFile.write('xticklabels at=edge bottom,xlabels at=edge bottom},');TeXFile.write('\n')
            TeXFile.write(r'ymajorgrids=true,xmajorgrids=true,ylabel=$\sigma_{12} \: (Pa)$,');TeXFile.write('\n')
            TeXFile.write('axis on top,scale only axis,width=0.4\linewidth,ymin=0,ymax='+str(ylim));TeXFile.write('\n')
            #TeXFile.write(', every x tick scale label/.style={at={(xticklabel* cs:1.05,0.75cm)},anchor=near yticklabel}]');TeXFile.write('\n')
            ## Black to white color map
            #TeXFile.write(', every x tick scale label/.style={at={(xticklabel* cs:1.05,0.75cm)},anchor=near yticklabel},colormap={bw}{gray(0cm)=(1); gray(1cm)=(0.05)}]');TeXFile.write('\n')
            ## Red to yellow color map
            TeXFile.write(', every x tick scale label/.style={at={(xticklabel* cs:1.05,0.75cm)},anchor=near yticklabel},colormap={ry}{rgb255(0cm)=(255,255,0);rgb255(1cm)=(255,0,0)}]');TeXFile.write('\n')
            ## Green to yellow color map
            #TeXFile.write(', every x tick scale label/.style={at={(xticklabel* cs:1.05,0.75cm)},anchor=near yticklabel},colormap={gy}{rgb255(0cm)=(255,255,0);rgb255(1cm)=(0,128,0)}]');TeXFile.write('\n')
            for k,rando in enumerate(nom):
                if k==0 :
                    TeXFile.write(r'\nextgroupplot[xlabel='+str(xlabels[i][k])+']')
                elif k==1 :
                    if rando[2:6]=='slow':
		        TeXFile.write(r'\nextgroupplot[colorbar,colorbar style={title= {$ c_s \: (m/s)$},every y tick scale label/.style={at={(2.,-.1125)}} },xlabel='+str(xlabels[i][k])+']')
                    elif rando[2:6]=='fast':
		        TeXFile.write(r'\nextgroupplot[colorbar,colorbar style={title= {$c_f \: (m/s)$},every y tick scale label/.style={at={(2.,-.1125)}} },xlabel='+str(xlabels[i][k])+']')
                TeXFile.write('\n')
                for j,name in enumerate(pgfFiles[i][k]):
                    
                    if name[31:31+6]=='_yield': ##  yield surface
                        TeXFile.write(r'\addplot[gray,dashed,thin] table[x='+str(srcX[k])+',y='+str(srcY[k])+'] {chapter5/pgfFigures/'+name+'};')
                        TeXFile.write('\n')
                    else:
                        TeXFile.write(r'\addplot[mesh,point meta = \thisrow{p},very thick,no markers] table[x='+str(srcX[k])+',y='+str(srcY[k])+'] {chapter5/pgfFigures/'+name+r'} node[above right,black] {$\textbf{'+str(j+1)+'}$};')
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
