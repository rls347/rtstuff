import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt

#Get the Henyey Greenstein coefficients chi_L
def hgcoef(n,g):
    hgout = np.zeros(n)
    for L in range(n):
        hgout[L] = (2*L+1)*(g**L)
    return hgout
    
def Phase(theta,g):
    pval = (1-g**2)/((1+(g**2)-(2*g*np.cos(theta)))**1.5)
    return pval

mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['font.size'] = 14

g=0.85
hg85 = hgcoef(181,g)
hg1 = hgcoef(100,1)
hg80 = hgcoef(100,0.8)
plt.plot(hg85,label='g = .85')
plt.plot(hg80,label='g = .8')
plt.legend()
plt.yscale('log')
plt.ylim(.001,100)
plt.title(r'Henyey Greenstein $\chi_\ell$')
plt.ylabel(r'$\chi_\ell$')
plt.xlabel(r'Order ($\ell$)')
plt.tight_layout()
plt.savefig('hgcoef.png')
plt.close()


angles = np.arange(181)  #radians
p = Phase(angles*np.pi/180.,g)
plt.plot(angles,p)
plt.plot(angles,hg85)
plt.yscale('log')
plt.ylim(.001,100)
plt.savefig('Phasemaybe.png')
plt.close()







