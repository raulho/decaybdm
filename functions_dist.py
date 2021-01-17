import numpy as np

t0=2.4e+19 # s

def time(z):
    return (2.4e+19)*np.power(1.0+z,-2.0)  # 1/s

def Hz(z):
    return (2.1e-20)*np.power(1.0+z,2.0)  # 1/s
    
def rhog(z):
    return 0.26*np.power(1.0+z,4.0)  # eV/cm3

def gamma(zx):
    return np.power(zx,2.0)/t0     # 1/s


#-------------------------------#
def rhoc(z,zx,vc):
    g_c=1.0/np.sqrt(1.0+np.power(vc,2.0))
    fa= np.sqrt(1.0+np.power(g_c,2.0)*np.power(vc,2.0)*np.power((1.0+z)/(1.0+zx),2.0))/g_c
    f0= np.sqrt(1.0+np.power(g_c,2.0)*np.power(vc,2.0)*np.power(1.0/(1.0+zx),2.0))/g_c
    return (3.9e+03)*np.power(1.0+z,3.0)*(fa/f0)  # eV/cm3

def rhoc1(z,zx,vc):
    return rhoc2(z,zx,vc)*np.power((1.0+z)/(1.0+zx),4.0)

def rhoc2(z,zx,vc):
    g_c=1.0/np.sqrt(1.0+np.power(vc,2.0))
    fa= np.sqrt(1.0+np.power(g_c,2.0)*np.power(vc,2.0)*np.power((1.0+z)/(1.0+zx),2.0))/g_c
    f0= np.sqrt(1.0+np.power(g_c,2.0)*np.power(vc,2.0)*np.power(1.0/(1.0+zx),2.0))/g_c
    return (3.9e+03)*np.power(1.0+z,3.0)*(fa/f0)  # eV/cm3
    
def Lmu(z):
    return (1.0-np.exp(-(np.power((1.0+z)/(5.8e+04),1.88))))*np.exp(-(np.power(z/2e+06,2.5)))

def bdm_dn1(vc,z,fx,zx):
    return rhoc1(z,zx,vc)*(fx*gamma(zx)*np.exp(-gamma(zx)*time(z)))/(Hz(z)*rhog(z)*(1.0+z))

def bdm_dn2(vc,z,fx,zx):
    return rhoc2(z,zx,vc)*(fx*gamma(zx)*np.exp(-gamma(zx)*time(z)))/(Hz(z)*rhog(z)*(1.0+z))

def func_intn1(vc,z,fx,zx):
    return bdm_dn1(vc,z,fx,zx)*Lmu(z)

def func_intn2(vc,z,fx,zx):
    return bdm_dn2(vc,z,fx,zx)*Lmu(z)

#-------------------------------#



