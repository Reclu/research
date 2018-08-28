#!/usr/bin/python

import numpy as np

def computeAnalyticalSolution(x,c,cp,t,tu,sigd,Sigy,HT,E):
    Sexact = np.zeros(len(x))
    EPexact = np.zeros(len(x))
    #The analytical solution is defined according to the current time t
    t1 = ((c*tu)/(c-cp))
    beta = (c+cp)/(c-cp)
    l1 = t1*cp
    t2 = t1 + l1/c
    t3 = (c*t2)/(c-cp)
    t4 = t3 + t3-t2
    S6 = 0.5*(sigd -Sigy)*((c/cp)-1.0) - (sigd/beta)
    EP1 = (sigd-Sigy)*((1.0/HT)-(1.0/E))
    if (t<t1):
        for i,valx in enumerate(x):
            if ((valx>(c*(t-tu))) and (valx<=(cp*t))):
                Sexact[i] = sigd
            elif ((valx>(cp*t)) and (valx<=(c*t))):
                Sexact[i] = Sigy
            if (valx<=(cp*t)):
                EPexact[i] = EP1
    if (t>t1) and (t<t2):
        for i,valx in enumerate(x):
            if ((valx>(c*(t2-t))) and (valx<=(cp*t))):
                Sexact[i] = sigd/beta
            elif ((valx>(cp*t)) and (valx<=(c*t))):
                Sexact[i] = Sigy
            if (valx<=(cp*t1)):
                EPexact[i] = EP1
            elif ((valx>(cp*t1)) and (valx<=cp*t)):
                EPexact[i] = ((sigd/beta)-Sigy)*((1.0/HT)-(1.0/E))
    if ((t>t2) and (t<t3)):
        for i,valx in enumerate(x):
            if ((valx>(c*(t-t2))) and (valx<=(cp*t))):
                Sexact[i] = sigd/beta
            elif ((valx>(cp*t)) and (valx<=(c*t))):
                Sexact[i] = Sigy
            if (valx<=(cp*t1)):
                EPexact[i] = EP1
            elif ((valx>(cp*t1)) and (valx<=(cp*t))):
                EPexact[i] = ((sigd/beta)-Sigy)*((1.0/HT)-(1.0/E))
    if ((t>t3) and (t<t4)):
        for i,valx in enumerate(x):
            if ((valx>(c*(t4-t))) and (valx<=(c*(t-t2)))):
                Sexact[i] = S6
            elif ((valx>(c*(t-t2))) and (valx<=(c*t))):
                Sexact[i] = Sigy
            if (valx<=(cp*t1)):
                EPexact[i] = EP1
            elif ((valx>(cp*t1)) and (valx<=cp*t3)):
                EPexact[i] = ((sigd/beta)-Sigy)*((1.0/HT)-(1.0/E))
    elif (t>t4):
        for i,valx in enumerate(x):
            if ((valx>(c*(t-t4))) and (valx<=(c*(t-t2)))):
                Sexact[i] = S6
            elif ((valx>(c*(t-t2))) and (valx<=(c*t))):
                Sexact[i] = Sigy
            if (valx<=(cp*t1)):
                EPexact[i] = EP1
            elif ((valx>(cp*t1)) and (valx<=(cp*t3))):
                EPexact[i] = ((sigd/beta)-Sigy)*((1.0/HT)-(1.0/E))
    return Sexact,EPexact
