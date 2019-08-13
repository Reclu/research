# !/usr/bin/python

import numpy as np

# Mesh definition
dx=1.
dy=1.

# Number of material point in the cell
N=4

# Material points positions
if N==1:
    xk=0.;yk=0.
elif N==4:
    xk=np.array([-0.25,+0.25,-0.25,0.25])
    yk=np.array([-0.25,-0.25,+0.25,+0.25])

# DCU method -> 3 cells 
S = np.zeros((12,3*N))
dSx = np.zeros((12,3*N))
dSy = np.zeros((12,3*N))

# Convention : Left element has number 1 // Current element has number 2 // Bottom element has number 3 // Bottom Left element has number 4
## Points are assumed to be at the same position in each cell
S1=(1-xk)*(1-yk)/4.          
S2=(1+xk)*(1-yk)/4.
S3=(1+xk)*(1+yk)/4.
S4=(1-xk)*(1+yk)/4.

S[0,0:N]=S1 ; S[1,0:N]=S2; S[2,0:N]=S3 ; S[3,0:N]=S4;
S[4,N:2*N]=S1 ; S[5,N:2*N]=S2; S[6,N:2*N]=S3 ; S[7,N:2*N]=S4;
S[8,2*N:3*N]=S1 ; S[9,2*N:3*N]=S2; S[10,2*N:3*N]=S3 ; S[11,2*N:3*N]=S4;

dS1xi=-(1-yk)/4. 
dS2xi= (1-yk)/4. 
dS3xi= (1+yk)/4. 
dS4xi=-(1+yk)/4. 

dS1et=-(1-xk)/4. 
dS2et=-(1+xk)/4. 
dS3et= (1+xk)/4. 
dS4et= (1-xk)/4. 

dSx[0,0:N]=dS1xi ; dSx[1,0:N]=dS2xi; dSx[2,0:N]=dS3xi ; dSx[3,0:N]=dS4xi;
dSx[4,N:2*N]=dS1xi ; dSx[5,N:2*N]=dS2xi; dSx[6,N:2*N]=dS3xi ; dSx[7,N:2*N]=dS4xi;
dSx[8,2*N:3*N]=dS1xi ; dSx[9,2*N:3*N]=dS2xi; dSx[10,2*N:3*N]=dS3xi ; dSx[11,2*N:3*N]=dS4xi;

dSy[0,0:N]=dS1et ; dSy[1,0:N]=dS2et; dSy[2,0:N]=dS3et ; dSy[3,0:N]=dS4et;
dSy[4,N:2*N]=dS1et ; dSy[5,N:2*N]=dS2et; dSy[6,N:2*N]=dS3et ; dSy[7,N:2*N]=dS4et;
dSy[8,2*N:3*N]=dS1et ; dSy[9,2*N:3*N]=dS2et; dSy[10,2*N:3*N]=dS3et ; dSy[11,2*N:3*N]=dS4et;

Sum=np.array([np.sum(S[0,:]),np.sum(S[1,:]),np.sum(S[2,:]),np.sum(S[3,:]),np.sum(S[0,:]),np.sum(S[1,:]),np.sum(S[2,:]),np.sum(S[3,:]),np.sum(S[0,:]),np.sum(S[1,:]),np.sum(S[2,:]),np.sum(S[3,:])]) ;

Phi1x=-(S[1,:]/Sum[1] + S[2,:]/Sum[2]) ;
Phi2x= (S[5,:]/Sum[5] + S[6,:]/Sum[6]) ;
Phi3x= (S[5,:]/Sum[5] + S[6,:]/Sum[6]) ;
Phi4x=-(S[1,:]/Sum[1] + S[2,:]/Sum[2]) ;

Phi1y=-(S[10,:]/Sum[10] + S[11,:]/Sum[11]) ;
Phi2y=-(S[10,:]/Sum[11] + S[11,:]/Sum[11]) ;
Phi3y= (S[6,:]/Sum[6] + S[7,:]/Sum[7]) ;
Phi4y= (S[6,:]/Sum[6] + S[7,:]/Sum[7]) ;

zero=np.zeros(3)

Phi_x = np.array([[np.zeros((4,3*N))],[Phi1x],[Phi2x],[Phi3x],[Phi4x],[np.zeros((4,3*N))]]);
Phi_y = np.array([np.zeros((4,3*N)),Phi1y,Phi2y,Phi3y,Phi4y,np.zeros((4,3*N))]);
print Phi_x.shape
"""
transverse = False

if transverse:
    points=4*N;
    Sum=[Sum,sum(S(1,:)),sum(S(2,:)),sum(S(3,:)),sum(S(4,:))] ;
    S=[S,zeros(size(S,1),N); zeros(4,4*N)];
    Phi_x = [Phi_x,zeros(size(Phi_x,1),N);zeros(4,4*N)];
    Phi_y = [Phi_y,zeros(size(Phi_y,1),N);zeros(4,4*N)];
    S(13,3*N+1:4*N)=S1 ; S(14,3*N+1:4*N)=S2; S(15,3*N+1:4*N)=S3 ; S(16,3*N+1:4*N)=S4;
    
    Phix_T=np.zeros(np.shape(Phi_x));
    Phiy_T=np.zeros(np.shape(Phi_y));
    
    Phix_T[4,:]=(S[0,:]/Sum[0]+S[1,:]/Sum[1]-S[14,:]/Sum[14]-S[15,:]/Sum[15])
    Phix_T[5,:]=-(S[4,:]/Sum[4]+S[5,:]/Sum[5]-S[10,:]/Sum[10]-S[11,:]/Sum[11])
    Phix_T[6,:]=-(S[4,:]/Sum[4]+S[5,:]/Sum[5]-S[10,:]/Sum[10]-S[11,:]/Sum[11])
    Phix_T[7,:]=(S[0,:]/Sum[0]+S[1,:]/Sum[1]-S[14,:]/Sum[14]-S[15,:]/Sum[15])
    
    Phiy_T[4,:]=(S[8,:]/Sum[8]+S[11,:]/Sum[11]-S[13,:]/Sum[13]-S[14,:]/Sum[14]) 
    Phiy_T[5,:]=(S[8,:]/Sum[8]+S[11,:]/Sum[11]-S[13,:]/Sum[13]-S[14,:]/Sum[14]) 
    Phiy_T[6,:]=-(S[4,:]/Sum[4]+S[7,:]/Sum[7]-S[1,:]/Sum[1]- S[2,:]/Sum[2]) 
    Phiy_T[7,:]=-(S[4,:]/Sum[4]+S[7,:]/Sum[7]-S[1,:]/Sum[1]- S[2,:]/Sum[2]) 
"""
