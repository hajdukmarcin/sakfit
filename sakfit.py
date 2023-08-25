import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
import matplotlib 
from scipy.optimize import curve_fit
import scipy.optimize as optimize
from numpy import sqrt, pi, exp, linspace

data = np.loadtxt('sakurai.dat')
x = data[:,0]
y = data[:,1]
ye = data[:,2]

a_min = 0.
a_max = 2000.
b_min = 0.1
b_max = 7.5
Steps = 251
Grid = np.zeros([Steps,Steps])
for s1 in range(Steps):
    for s2 in range(Steps):
        a = a_min + (a_max - a_min)*float(s1)/float(Steps-1)
        b = b_min + (b_max - b_min)*float(s2)/float(Steps-1)
        chi2 = 0.0
        for n in range(len(x)):
            residual = y[n] - 0.05436*a*2.77e-1*x[n]*x[n]*(1.-np.exp(-b*13.26*2.1e5*np.power(a*10000,-1.35)*(np.power(x[n],-2.1))))
            chi2 = chi2 + residual*residual/(ye[n]*ye[n])
        rchi2 = chi2 / (n - 1.0)
        Grid[Steps-1-s2,s1] = rchi2

plt.subplot(211)
plt.figure(1, figsize=(4,6))
mini  = np.min(Grid)  # minimal value of chi2
np.argmin(Grid)
image = plt.imshow(Grid, vmin=mini, vmax=mini+1.5, extent=[a_min*17.28,a_max*17.28,b_min*4003.26,b_max*4003.26])
plt.colorbar(image)

plt.ylabel(r'$EM [10^{5} cm^{-6} pc]$', fontsize=10)
plt.xlabel(r'$T_e [k]$', fontsize=10)
plt.show()
