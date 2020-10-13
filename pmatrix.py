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

vals = np.polynomial.legendre.leggauss(nangles)  
quad = vals[0]
weights = vals[1]

quadneg = quad[:8]
quadpos = quadneg*-1

chi = hgcoef(32,g)

pleg = np.zeros((32,8))
plegneg = np.zeros((32,8))

for N in range(32):
    pleg[N,:] = P(N,quadpos)
    plegneg[N,:] = P(N,quadneg)

print('positive pleg')
print(pleg)
print(' ')
print('negative pleg')
print(plegneg)
print(' ')
#This is what Graeme has written as 'solar phase vector' so this is just for each angle wrt the sun
print('Positive Phase', phasefunc(quadpos,32,g))
print('Negative Phase', phasefunc(quadneg,32,g))
print(' ')

posphasematrix = np.zeros((nquad,nquad))
negphasematrix = np.zeros((nquad,nquad))

#pos phase matrix is plus angles with plus angles or neg with neg
#neg is plus with neg
for i in range(8):
    for j in range(8):
        posphasematrix[i,j] = scatterphase(pleg,pleg,i,j,32,g)
        negphasematrix[i,j] = scatterphase(pleg,plegneg,i,j,32,g)

print('P+ phase matrix')
print(posphasematrix)
print(' ')
print('P- phase matrix')
print(negphasematrix)

#rmatrix = (1+delta(0,m)) * (omega/4) * M**(-1m) * P- *C
#we still are doing only m = 0, so delta = 1, the M part is just M
#M is diagonal with mus, C is diagonal with weights

#const = omega/2.
#what is omega? single scattering albedo... ? 

M = np.diag(quad[:8])
C = np.diag(weights[:8])

print('r matrix')
rmatrix = np.matrix(M)*np.matrix(negphasematrix)*np.matrix(C)
print(rmatrix)

print(M)
print(C)


