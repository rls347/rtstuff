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

g=.85
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
plegsolar=np.zeros((32,1))

for N in range(32):
    pleg[N,:] = P(N,quadpos)
    plegneg[N,:] = P(N,quadneg)
    plegsolar[N,0] = P(N,solar)

print('positive pleg')
print(pleg)
print(' ')
print('negative pleg')
print(plegneg)
print(' ')
print('solar pleg')
print(plegsolar)
print(' ')

possun = np.zeros((nquad,1))
negsun = np.zeros((nquad,1))

for i in range(8):
    for j in range(1):
        possun[i,j] = scatterphase(pleg,plegsolar,i,j,32,g)
        negsun[i,j] = scatterphase(plegneg,plegsolar,i,j,32,g)

print('Positive Solar Phase Vector')
print(possun)
print('Negative Solar Phase Vector')
print(negsun)
print(' ')

posphasematrix = np.zeros((nquad,nquad))
negphasematrix = np.zeros((nquad,nquad))

#pos phase matrix is plus angles with plus angles or neg with neg
#neg is plus with neg
for i in range(8):
    for j in range(8):
        posphasematrix[i,j] = scatterphase(pleg,pleg,i,j,32,g)
        negphasematrix[i,j] = scatterphase(plegneg,pleg,i,j,32,g)

print('P+ phase matrix')
print(posphasematrix)
print(' ')
print('P- phase matrix')
print(negphasematrix)

#rmatrix = (1+delta(0,m)) * (omega/4) * M**(-1m) * P- *C
#we still are doing only m = 0, so delta = 1, the M part is just M
#M is diagonal with mus, C is diagonal with weights

#const = omega/2.
#what is omega? single scattering albedo... ? needs to be 1 to get right r matrix

mm = 1./(quad[8:])
cc = weights[8:]

M = np.matrix(np.diag(mm[::-1]))
C = np.matrix(np.diag(cc[::-1]))

print('r matrix')
rmatrix = .5*M*np.matrix(negphasematrix)*C
print(rmatrix)

# print(M)
# print(C)
#
# print(quad)
# print(weights)

print('')
print('t matrix')
tmatrix = M- .5*M*np.matrix(posphasematrix)*C
print(tmatrix)

#plus/minus versus minus/plus in Graeme's notes. 
#F/4pi e(-tau/mu) don't know what F should be but as 1 it works out I think
#Yes - Graeme has Ftop = 1.0
#in his solution tau = 1. not sure why optical depth 1 makes any sense here. 
#mu=cos(0)=1 mu is mu(dot) here. sun angle. 

const = np.exp(-1.)/(4.*np.pi)

sourcepos = const * M * np.matrix(negsun) #- (1-omega0)*Beta*M*U
sourceneg = const * M * np.matrix(possun)

print('source pos')
print(sourcepos)
print('source neg')
print(sourceneg)
print(' ')










