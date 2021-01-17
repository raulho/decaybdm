#--------------------------------------------------------------
# Code for getting CMB spectral μ-distortions on
# BDM model modified
# Autor: Raúl Henríquez
# date: 21/09/2020
#--------------------------------------------------------------

import numpy as np
import scipy.integrate as integrate
from functions_dist import *


N=50

ac=np.linspace(4.9e-07,3.36e-06,N)
vc=np.linspace(0,0.71,N)
zc=(1.0/ac)-1.0

feff=1.0
fx=1.0

i0=0.0

#-------------------------------------------#
f = open('data_mu_fx_feff.dat', 'w')

for i in range(N):
    for j in range(N):
        invexp1 = lambda x: func_intn2(vc[j],x,fx,zc[i])
        y1,err=integrate.quad(invexp1,0.0,zc[i])       

        invexp2 = lambda x: func_intn1(vc[j],x,fx,zc[i])
        y2,err=integrate.quad(invexp2,zc[i],np.inf)      
    
        mu_bdm = 1.401*feff*(y1+y2*i0)

        f.write("%8f %.9f %.12f\n" % (zc[i], vc[j], mu_bdm))


f.close()
#-------------------------------------------#              

 
