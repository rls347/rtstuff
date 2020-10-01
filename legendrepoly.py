import numpy as np
import matplotlib.pyplot as plt

def P(n,x):
    if n==0:
        return 1
    elif n==1:
        return x
    elif n>1:
        return (1./n)* ((2*(n-1)+1)*x*P(n-1,x) - (n-1)*P(n-2,x))
        
        
def testP6(x):
    return 1./16 * (231*x**6 - 315*x**4 + 105*x**2 -5)

def testP4(x):
    return 1./8 * ( 35*x**4 - 30*x**2 + 3)
    
def testP2(x):
    return 1./2 * (3*x**2 - 1)

# PN(1) = 1
# P0(x) = 1
# P1(x) = x
# (n+1)*P[n+1] = (2*n+1)*x*P[n] - n*P[n-1]
# (n)*P[n] = (2*n)*x*P[n-1] - (n-1)*P[n-2]

# Test - from Wolfram
# P6(x) = 1/16 * (231x**6 - 315x**4 + 105x**2 -5)
# PN(1) = 1. easy test
# print(testP6(1),P(6,1))
# print(testP4(1),P(4,1))
# print(testP2(1),P(2,1))
#
# for i in range(7):
#     print(i,P(i,1))

angles = np.arange(181)  #radians
mu = np.cos(angles*np.pi/180.)

plt.plot(mu,P(7,mu),label='N=7')
plt.plot(mu,P(15,mu),label='N=15')
plt.plot(mu,P(23,mu),label='N=23')
plt.plot(mu,P(31,mu),label='N=31')
plt.legend()
plt.savefig('polynomials.png')
plt.close()













