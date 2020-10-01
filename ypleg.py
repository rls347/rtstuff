import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt

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

mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['font.size'] = 14

def P(n,x):
    if n==0:
        return 1
    elif n==1:
        return x
    elif n>1:
        return (1./n)* ((2*(n-1)+1)*x*P(n-1,x) - (n-1)*P(n-2,x))
 
#testing graeme's version. moving his code over but we only use M=0 (??)
def ypleg(n,mu):
    if n==0:
        return 1
    elif n==1:
        return mu
    else:
        L=n-2
        LM2 = (L+2)*(L+2)
        LM1 = (L+1)*(L+1)
        #print(n,L,LM2,LM1)
        return ((2*(n-1)+1) * mu * ypleg(n-1,mu) - np.sqrt(LM1)*ypleg(n-2,mu))/np.sqrt(LM2)
        #return N=L+2, L=N-2
        #YPLEG(N+2,I) = ((2*(L+1)+1)*MU*YPLEG(L+1,I)-SQRT(LM1)*
        #     $            YPLEG(L,I))/SQRT(LM2)
    
    

angles = np.arange(181)  #radians
mu = np.cos(angles*np.pi/180.)

phases = np.zeros_like(mu)
n=7
g=.8

# p7 = phasefunc(mu,7,.8)
# p15 = phasefunc(mu,15,.8)
# p23 = phasefunc(mu,23,.8)
# p31 = phasefunc(mu,31,.8)
#
# plt.plot(angles,p7,label='n=7')
# plt.plot(angles,p15,label='n=15')
# plt.plot(angles,p23,label='n=23')
# plt.plot(angles,p31,label='n=31')
# #plt.plot(angles,Phase1(mu,g),label='P')
# plt.yscale('log')
# plt.ylim(.001,100)
# plt.title('Phase Functions')
# plt.xlabel('Scattering Angle')
# plt.legend()
# plt.tight_layout()
# plt.savefig('phasefunction.png')
# plt.close()


for n in range(20):
    print(ypleg(n,.5),P(n,.5))


    
