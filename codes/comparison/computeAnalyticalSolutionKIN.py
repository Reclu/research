#!/usr/bin/python

import numpy as np
import pdb

def computeAnalyticalSolutionKIN(x,L,c,t,vd,HEL,E,H,rho):
    #Definition and initialization of arrays and data
    lam=mu=1.
    Sexact = np.zeros(len(x))
    EPexact = np.zeros(len(x))
    Vexact = np.zeros(len(x)) 
    HT =  E*H/(E+H)#(2.0*lam*mu+(lam+2.0*mu)*(mu+3.0*H/2.0))/(3.0*mu+3.0*H/2.0)
    cp = np.sqrt(HT/rho)
    ##Comment: The analytical solution is defined according to the current time t
    #times (and positions) associated to wave interactions
    t1 = (L/2.0)/c
    t2 = L/(c+cp);x2p = L+(t2/2.0)*(cp-c);x2m = (t2/2.0)*(c-cp)
    t3 = (L-x2p)/c+t2
    t4 = (x2p-x2m)/(2.0*c) +t2
    t5 = L/(4.0*c) + (t3+t4)/2.0 ; x5p = 3.0*L/4.0 + c*(t3-t4)/2.0
    x5m = L/2.0- (x5p-L/2.0)
    t6 = (x5p-L/2.0)/c + t5
    t7 = (t5+t6)/2.0 + (x5p-(L/2.0))/(2.0*cp) 
    x7p = (x5p+(L/2.0))/2.0 + cp*(t5-t6)/2.0
    x7m = L/2.0- (x7p-L/2.0) 
    t8 = (L/2.0-x7m)/cp + t7
    tint = (L/2.0-x7m)/cp + t8
    t9 = (L-x5p)/c+t5
    #Stress and velocity states
    v1 = HEL/(rho*c)-vd ; v1p = -v1 ; v2 = 0.0
    v3 = 2.0*HEL/(rho*c) - vd ; v3p = -v3
    S2 =  HEL*(1.0-(cp/c))+rho*cp*vd 
    v4 = S2/(2.0*rho*c) + v3/2.0 ; v4p = -v4
    S4 = (S2-(rho*c*v3))/2.0
    S5 = -2.0*HEL+rho*c*vd
    S7 = HEL*(0.2*cp/c - 1.0)
    v4 = S2/(2.0*rho*c) + v3/2.0 ;  v6 =  S4/(rho*c) + v4 ; v6p = -v6 ; v5 = 0.0
    v7 = S7/(rho*c)+v6 ; v7p = (S5-S7)/(rho*c) ; v7s = -v7p ; v7t = -v7
    S8 = S7 + rho*cp*(v7p-v7)/2.0
    v8 = (v7+v7p)/2.0 ; v8p = -v8
    S9 = S7 - rho*cp*v7p ; v9 = (v7p+v7s)/2.0
    S10 = (S8+S9-rho*cp*v8)/2.0
    v10 = (S9-S8)/(2.0*rho*cp) + v8/2.0 ; v10p = -v10
    S11 = S10 - rho*cp*v10 ; v11 = 0.0
    #Plastic strain states
    EP2 = ((S2 - HEL)/(2.0*mu))*(((E)/HT)-1.0)
    EP2 = (S2-HEL)*((1.0/HT)-(1.0/E))
    KK = mu*(3.0*lam+2.0*mu)/(lam+2.0*mu) + 3.0*H/2.0
    EP8 = 2*mu*(S8+HEL)/(KK*(lam+2.0*mu))
    EP9 = 2*mu*(S9+HEL)/(KK*(lam+2.0*mu))
    EP10 = 2*mu*(S10+HEL)/(KK*(lam+2.0*mu))
    EP11 =2*mu*(S11+HEL)/(KK*(lam+2.0*mu))
    if (t<t1):
        for i,valx in enumerate(x):
            if ((valx>=(-c*t+(L/2.0))) and (valx<=(-cp*t+(L/2.0)))):
                Sexact[i] = HEL
                Vexact[i] = v1
            elif ((valx>=(cp*t+(L/2.0))) and (valx<=(c*t+(L/2.0)))):
                Sexact[i] = HEL
                Vexact[i] = v1p
            elif ((valx>=(-cp*t+(L/2.0))) and (valx<=(cp*t+(L/2.0)))):
                Sexact[i] = S2
                Vexact[i] = v2
                EPexact[i] = EP2
            elif (valx<(-c*t+(L/2.0))):
                Vexact[i] = -vd
            elif (valx>(c*t+(L/2.0))):
                Vexact[i] = vd
    elif (t>=t1) and (t<t2):
        for i,valx in enumerate(x):
            if ((valx>(c*(t-t1))) and (valx<(-cp*t+(L/2.0)))):
                Sexact[i] = HEL
                Vexact[i] = v1
            elif ((valx<(L-c*(t-t1))) and (valx>(cp*t+(L/2.0)))):
                Sexact[i] = HEL
                Vexact[i] = v1p
            elif ((valx>=(-cp*t+(L/2.0))) and (valx<=(cp*t+(L/2.0)))):
                Sexact[i] = S2
                Vexact[i] = v2
                EPexact[i] = EP2
            elif (valx<(c*(t-t1))):
                Vexact[i] = v3     
            elif (valx>(L-c*(t-t1))):
                Vexact[i] = v3p  
    elif ((t>=t2) and (t<=t3)):
        for i,valx in enumerate(x):
            if ((valx>x2m-c*(t-t2)) and (valx<x2m+c*(t-t2))):
                Sexact[i] = S4
                Vexact[i] = v4
            elif ((valx<x2p+c*(t-t2)) and (valx>x2p-c*(t-t2))):
                Sexact[i] = S4
                Vexact[i] = v4p
            elif (valx>x2m+c*(t-t2)) and (valx<x2p-c*(t-t2)):
                Sexact[i] = S2  
                Vexact[i] = v2
            elif (valx<x2m-c*(t-t2)):
                Vexact[i] = v3     
            elif (valx>x2p+c*(t-t2)):
                Vexact[i] = v3p  
            if ((valx>=x2m) and (valx<=x2p)):
                EPexact[i] = EP2
    elif ((t>t3) and (t<=t4)):
        for i,valx in enumerate(x):
            if ((valx>c*(t-t3)) and (valx<x2m+c*(t-t2))):
                Sexact[i] = S4 
                Vexact[i] = v4
            elif ((valx<L-c*(t-t3)) and (valx>x2p-c*(t-t2))):
                Sexact[i] = S4 
                Vexact[i] = v4p
            elif (valx>x2m+c*(t-t2)) and (valx<x2p-c*(t-t2)):
                Sexact[i] = S2 
                Vexact[i] = v2
            elif (valx<c*(t-t3)):
                Vexact[i] = v6 
            elif (valx>L-c*(t-t3)):
                Vexact[i] = v6p       
            if ((valx>=x2m) and (valx<=x2p)):
                EPexact[i] = EP2
    elif ((t>t4) and (t<t5)):
        for i,valx in enumerate(x):
            if ((valx>c*(t-t3)) and (valx<L/2.0-c*(t-t4))):
                Sexact[i] = S4   
                Vexact[i] = v4
            elif ((valx<L-c*(t-t3)) and (valx>L/2.0+c*(t-t4))):
                Sexact[i] = S4 
                Vexact[i] = v4p
            elif (valx>L/2.0-c*(t-t4)) and (valx<L/2.0+c*(t-t4)):
                Sexact[i] = S5
                Vexact[i] = v5
            elif (valx<c*(t-t3)):
                Vexact[i] = v6 
            elif (valx>L-c*(t-t3)):
                Vexact[i] = v6p 
            if ((valx>=x2m) and (valx<=x2p)):
                EPexact[i] = EP2
    elif ((t>=t5) and (t<t6)):
        for i,valx in enumerate(x):
            if ((valx>x5m-c*(t-t5)) and (valx<x5m-cp*(t-t5))):
                Sexact[i] = S7
                Vexact[i] = v7
            elif ((valx>x5m+cp*(t-t5)) and (valx<x5m+c*(t-t5))):
                Sexact[i] = S7
                Vexact[i] = v7p
            elif ((valx>x5p-c*(t-t5)) and (valx<x5p-cp*(t-t5))):
                Sexact[i] = S7
                Vexact[i] = v7s
            elif ((valx>x5p+cp*(t-t5)) and (valx<x5p+c*(t-t5))):
                Sexact[i] = S7
                Vexact[i] = v7t
            elif ((valx>x5m-cp*(t-t5)) and (valx<x5m+cp*(t-t5))):
                Sexact[i] = S8
                Vexact[i] = v8
            elif ((valx>x5p-cp*(t-t5)) and (valx<x5p+cp*(t-t5))):
                Sexact[i] = S8
                Vexact[i] = v8p
            elif ((valx>x5m+c*(t-t5)) and (valx<x5p-c*(t-t5))):
                Sexact[i] = S5 
                Vexact[i] = v5
            elif (valx<x5m-c*(t-t5)):
                Vexact[i] = v6 
            elif (valx>x5p+c*(t-t5)):
                Vexact[i] = v6p 
            if (((valx>x5m-cp*(t-t5)) and (valx<x5m+cp*(t-t5))) \
                or ((valx>x5p-cp*(t-t5)) and (valx<x5p+cp*(t-t5)))):
                EPexact[i] = EP8
            elif (((valx>=x2m) and (valx<=x5m-cp*(t-t5))) \
               or ((valx>x5m+cp*(t-t5)) and (valx<x5p-cp*(t-t5))) \
               or  ((valx>x5p+cp*(t-t5)) and (valx<=x2p))):
                EPexact[i] = EP2
    elif ((t>=t6) and (t<t7)):
        for i,valx in enumerate(x):
            if ((valx>x5m-c*(t-t5)) and (valx<x5m-cp*(t-t5))):
                Sexact[i] = S7
                Vexact[i] = v7
            elif ((valx>x5m+cp*(t-t5)) and (valx<L/2.0-cp*(t-t6))):
                Sexact[i] = S7
                Vexact[i] = v7p
            elif ((valx>=L/2.0-cp*(t-t6)) and (valx<=L/2.0+cp*(t-t6))):
                Sexact[i] = S9
                Vexact[i] = v9           
            elif ((valx>L/2.0+cp*(t-t6)) and (valx<x5p-cp*(t-t5))):
                Sexact[i] = S7
                Vexact[i] = v7s
            elif ((valx>x5p+cp*(t-t5)) and (valx<x5p+c*(t-t5))):
                Sexact[i] = S7
                Vexact[i] = v7t
            elif ((valx>x5m-cp*(t-t5)) and (valx<x5m+cp*(t-t5))):
                Sexact[i] = S8
                Vexact[i] = v8
            elif ((valx>x5p-cp*(t-t5)) and (valx<x5p+cp*(t-t5))):
                Sexact[i] = S8
                Vexact[i] = v8p
            elif (valx<x5m-c*(t-t5)):
                Vexact[i] = v6 
            elif (valx>x5p+c*(t-t5)):
                Vexact[i] = v6p
            if (((valx>x5m-cp*(t-t5)) and (valx<x5m+cp*(t-t5))) \
                or ((valx>x5p-cp*(t-t5)) and (valx<x5p+cp*(t-t5)))):
                EPexact[i] = EP8
            elif ((valx>=L/2.0-cp*(t-t6)) and (valx<=L/2.0+cp*(t-t6))):
                EPexact[i] = EP9      
            elif (((valx>=x2m) and (valx<=x5m-cp*(t-t5))) \
                  or ((valx>=x5m+cp*(t-t5)) and (valx<L/2.0-cp*(t-t6))) \
                  or ((valx>L/2.0+cp*(t-t6)) and  (valx<=x5p-cp*(t-t5))) \
               or  ((valx>x5p+cp*(t-t5)) and (valx<=x2p))):
                EPexact[i] = EP2
    elif ((t>=t7) and (t<t8)):
        for i,valx in enumerate(x):
            if ((valx>x5m-c*(t-t5)) and (valx<x5m-cp*(t-t5))):
                Sexact[i] = S7
                Vexact[i] = v7
            elif ((valx>x5p+cp*(t-t5)) and (valx<x5p+c*(t-t5))):
                Sexact[i] = S7
                Vexact[i] = v7t
            elif ((valx>=x7m+cp*(t-t7)) and (valx<x7p-cp*(t-t7))):
                Sexact[i] = S9
                Vexact[i] = v9 
            elif ( valx>x7m-cp*(t-t7)) and (valx<x7m+cp*(t-t7)):
                Sexact[i] = S10
                Vexact[i] = v10
            elif ( valx>x7p-cp*(t-t7)) and (valx<x7p+cp*(t-t7)):
                Sexact[i] = S10
                Vexact[i] = v10p
            elif ((valx>x5m-cp*(t-t5)) and (valx<x7m-cp*(t-t7))):
                Sexact[i] = S8
                Vexact[i] = v8
            elif ((valx>x7p+cp*(t-t7)) and (valx<x5p+cp*(t-t5))):
                Sexact[i] = S8
                Vexact[i] = v8p
            elif (valx<x5m-c*(t-t5)):
                Vexact[i] = v6 
            elif (valx>x5p+c*(t-t5)):
                Vexact[i] = v6p
            if (((valx>x5m-cp*(t-t5)) and (valx<x7m-cp*(t-t7))) \
                or ((valx>x7p+cp*(t-t7)) and (valx<x5p+cp*(t-t5)))):
                EPexact[i] = EP8
            elif (((valx >x2m) and (valx<=x5m-cp*(t-t5))) \
               or ((valx>=x5p+cp*(t-t5)) and (valx<x2p))):
                EPexact[i] = EP2
            elif ((valx>x7m+cp*(t-t7)) and (valx<x7p-cp*(t-t7))):
                EPexact[i] = EP9      
            elif ((( valx>x7m-cp*(t-t7)) and (valx<=x7m)) \
              or (( valx>=x7p) and (valx<x7p+cp*(t-t7)))):
                EPexact[i] = EP10
            elif ((( valx>x7m) and (valx<x7m+cp*(t-t7))) \
                  or (( valx>x7p-cp*(t-t7)) and (valx<=x7p))):
                EPexact[i] = EP10
    elif ((t>=t8) and (t<t9)):
        for i,valx in enumerate(x):
            if ((valx>x5m-c*(t-t5)) and (valx<x5m-cp*(t-t5))):
                Sexact[i] = S7
                Vexact[i] = v7
            elif ((valx>x5p+cp*(t-t5)) and (valx<x5p+c*(t-t5))):
                Sexact[i] = S7
                Vexact[i] = v7t
            elif ((valx>x5m-cp*(t-t5)) and (valx<x7m-cp*(t-t7))):
                Sexact[i] = S8
                Vexact[i] = v8
            elif ((valx>x7p+cp*(t-t7)) and (valx<x5p+cp*(t-t5))):
                Sexact[i] = S8
                Vexact[i] = v8p
            elif (valx<x5m-c*(t-t5)):
                Vexact[i] = v6 
            elif (valx>x5p+c*(t-t5)):
                Vexact[i] = v6p
            elif ( valx>x7m-cp*(t-t7)) and (valx<L/2.0-cp*(t-t8)):
                Sexact[i] = S10
                Vexact[i] = v10
            elif ( valx>L/2.0+cp*(t-t8)) and (valx<x7p+cp*(t-t7)):
                Sexact[i] = S10
                Vexact[i] = v10p
            elif (valx>=L/2.0-cp*(t-t8)) and (valx<=L/2.0+cp*(t-t8)):
                Sexact[i] = S11
                Vexact[i] = v11  
            #
            if (((valx>x5m-cp*(t-t5)) and (valx<x7m-cp*(t-t7))) \
                or ((valx>x7p+cp*(t-t7)) and (valx<x5p+cp*(t-t5)))):
                EPexact[i] = EP8
            elif (((valx >x2m) and (valx<=x5m-cp*(t-t5))) \
               or ((valx>=x5p+cp*(t-t5)) and (valx<x2p))):
                EPexact[i] = EP2
            elif ((valx>=L/2.0-cp*(t-t8)) and (valx<=L/2.0+cp*(t-t8))):
                EPexact[i] = EP11      
            elif (( valx>x7m-cp*(t-t7) and (valx<L/2.0-cp*(t-t8))) \
                or (((valx>L/2.0+cp*(t-t8))) and  (valx<x7p+cp*(t-t7)))):
                EPexact[i] = EP10
    return Sexact,EPexact,Vexact

