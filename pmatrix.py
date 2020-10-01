import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import os

def hgcoef(n,g):
    hgout = np.zeros(n)
    for L in range(n):
        hgout[L] = (2*L+1)*(g**L)
    return hgout
    
def Phase1(mu,g):
    pval = (1-g**2)/((1+(g**2)-(2*g*mu))**1.5)
    return pval
    
def phasefunc(mu,n,g):
    phases = np.zeros_like(mu)
    chi = hgcoef(n,g)
    for l in range(n):
        phases= phases + chi[l]*P(l,mu)#*Phase1(mu,g)
    return phases    

def P(n,x):
    if n==0:
        if type(x) is np.ndarray:
            return(np.ones_like(x))
        else:
            return 1
    elif n==1:
        return x
    elif n>1:
        return (1./n)* ((2*(n-1)+1)*x*P(n-1,x) - (n-1)*P(n-2,x))
        
def scatterphase(p1,p2,i,j,n,g):
    chi = hgcoef(n,g)
    phase = 0.
    for l in range(n):
        phase = phase + chi[l]*p1[l,i]*p2[l,j]
    return phase
    
    
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['font.size'] = 14
#np.set_printoptions(linewidth=os.get_terminal_size().columns)
np.set_printoptions(linewidth=200)

g=.84
solar = 0.9 #cos of solar zenith angle

nquad=8
nangles=2*nquad
quad =np.zeros(nangles)
for i in range(1,nangles+1):
    quad[i-1] = np.cos(np.pi*(i-.25)/(nangles+.5))

#values from Graeme's text file. off only by 10-4    
Gr=np.array([-0.98940, -0.94458, -0.86563, -0.75540, -0.61788, -0.45802, -0.28160, -0.09501,
             0.09501, 0.28160, 0.45802, 0.61788, 0.75540, 0.86563, 0.94458, 0.98940])
quad = Gr
quadneg = Gr[:8]
quadpos = quadneg*-1

#weights = phasefunc(quad,n,g)
chi = hgcoef(32,g)

# pleg = np.zeros((32,8))
# plegneg = np.zeros((32,8))
#
# for N in range(32):
#     pleg[N,:] = P(N,quadpos)
#     plegneg[N,:] = P(N,quadneg)
#
# print('positive pleg')
# print(pleg)
# print(' ')
# print('negative pleg')
# print(plegneg)
# print(' ')
# #This is what Graeme has written as 'solar phase vector' so this is just for each angle wrt the sun
# print('Positive Phase', phasefunc(quadpos,32,g))
# print('Negative Phase', phasefunc(quadneg,32,g))
# print(' ')
#
# posphasematrix = np.zeros((nquad,nquad))
# negphasematrix = np.zeros((nquad,nquad))
# #
# # negneg = np.zeros((nquad,nquad))
# # negpos = np.zeros((nquad,nquad))
#
# #pos phase matrix is plus angles with plus angles or neg with neg
# #neg is plus with neg
# for i in range(8):
#     for j in range(8):
#         posphasematrix[i,j] = scatterphase(pleg,pleg,i,j,32,g)
#         negphasematrix[i,j] = scatterphase(pleg,plegneg,i,j,32,g)
#         # negneg[i,j] = scatterphase(plegneg,plegneg,i,j,32,g)
#         # negpos[i,j] = scatterphase(plegneg,pleg,i,j,32,g)
# print('P+ phase matrix')
# print(posphasematrix)
# print(' ')
# print('P- phase matrix')
# print(negphasematrix)

eps=3.0e-14

Z=quadpos[0]
for i in range(1,9):
    Z=quadpos[i-1]
    print(Z)
    for k in range(10):
        #print(k)
        p1=1.
        p2=0.
        for j in range(1,32):
            p3=p2
            p2=p1
            p1 = ((2*j-1)*Z*p2-(j-1*p3))/j
            print(p3,p2,p1)
        PP=16*(Z*p1-p2)/(Z*Z-1)
        Z1=Z
        Z=Z1-p1/PP
        print(k,PP)
        if np.abs(Z-Z1)>=eps:
            w=2./((1.-Z*Z)*PP*PP)
            k=10
    
#       N=2*NZ
#       DO 12 I=1,NZ
#         Z=DCOS(3.141592654D0*(I-.25D0)/(N+.5D0))
#         K=0
# 1       CONTINUE
#           P1=1.D0
#           P2=0.D0
#           DO 11 J=1,N
#             P3=P2
#             P2=P1
#             P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
# 11        CONTINUE
#           PP=N*(Z*P1-P2)/(Z*Z-1.D0)
#           Z1=Z
#           Z=Z1-P1/PP
#           K=K+1
#         IF (ABS(Z-Z1).GT.EPS .AND. K.LT.10) GO TO 1
#         ABSCISSAS(NZ+1-I)=Z
#         WEIGHTS(NZ+1-I)=2.D0/((1.D0-Z*Z)*PP*PP)
# 12    CONTINUE
#







