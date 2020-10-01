import numpy as np

#def PLEG(angleo,gmu,m,nexp,nzen,ypleg):

# g=0.85
# cosine of the solar zenith angle =0.9
# Gaussion quadrature (single); order 8+8
#
# S =            16
#  G =    0.8400000000000000
#  MU0 =     1.000000000000000
#  OMEGA =     1.000000000000000
#  FTOP =     1.000000000000000
#  TAU =     1.000000000000000
#
# nzen=8
# nexp=8
# angle0 = .9
# gmu = np.ones(nzen)
# ypleg = np.zeros((nexp,nzen+2))
# something=8
# for m in range(something):
#     if m==0:
#         cm1=1.
#         dm1=1.
#     else:
#         cm2 = cm1 * ((2*m-1)/(2*m))**(.5)
#         dm2 = dm1 * ((2*m+1)/(2*m))**(.5)
#
#     for i in range(nzen+2):
#         if i == nzen+1:
#             mu = 0
#         elif i==nzen:
#             mu = angle0
#         else:
#             mu = gmu[i-1]
#
#         if m==0:
#             ypleg[0,i] = 1.
#             ypleg[1,i] = mu
#         else:
#             ypleg[m,i] = cm2 * (1-mu**2)**(m/2.)
#             ypleg[m+1,i] = dm2 * mu * (1-mu**2)**(m/2.)
#
#         for il in range(1,nexp-1):
#             l=il-1
#             if l < m:
#                 ypleg[il,i] = 0
#             else:
#                 lm2 = (l-m+2)*(l+m+2)
#                 lm1 = (l-m+1)*(l+m+1)
#                 ypleg[il+1,i] = (((2*l+3)*mu*ypleg[il,i]) - ((lm1**.5) * ypleg[il-1,i])) / (lm2**.5)
#     if m>0:
#         cm1=cm2
#         dm1=dm2
#
# print(ypleg)
#
#
# def P(g,mu):
#     return (1-g**2)/( (1+g**2-2*g*mu)**(3./2))
#  N=2*NZ
#       DO 12 I=1,NZ
#         Z=DCOS(3.141592654D0*(I-.25D0)/(N+.5D0))
#
