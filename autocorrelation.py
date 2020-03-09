import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pickle
import numpy as np
import scipy.special as sp
import sys

#LaTeX
plt.rc('text', usetex=True)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 20}
plt.rc('font', **font)

# Import the parameter from param.py
import param

Nt_list = param.Nt_list
default_beta = param.default_beta
dr0 = param.dr0
N = param.N
precision = param.precision
f = param.f
choosen_beta = param.choosen_beta
real = param.real
imag = param.imag
beta_min = param.beta_min
beta_max = param.beta_max
beta_nb  = param.beta_nb
marker_list = param.marker_list
color_list = param.color_list

beta = np.linspace(beta_min,beta_max,beta_nb)

# Initialize list before loading data
expect = []
expect_with = []

expect_full = []
expect_with_full = []
######################## Function stat #################################

def Moyenne(ensemble):
    result = 0
    N = len(ensemble)
    for el in ensemble:
        result += el
    return result/N

def Variance(ensemble):
    result = 0
    N = len(ensemble)
    moy = Moyenne(ensemble)
    for el in ensemble:
        result += (el - moy)**2
    return result * (N-1)/N
######################### Reading the data #############################
for i in range(1000):
    with open("data/expectation_{0}_False.txt".format(i),"rb") as fp:
        expect_full.append(pickle.load(fp))
    with open("data/expectation_{0}_True.txt".format(i),"rb") as fp:
        expect_with_full.append(pickle.load(fp))

beta = 1.0
Nt = 20
N=len(expect_full)

sigma = []
sigma_with = []

x = np.linspace(1,N,N)

for Nbins in x:

    a=0
    for i in range(len(expect_full)):
        if i%Nbins!=0 or i==0:
            a+=expect_full[i]
        else:
            a = a/Nbins
            expect.append(a)
            a = 0

    a=0
    for i in range(len(expect_with_full)):
        if i%Nbins!=0 or i==0:
            a+=expect_with_full[i]
        else:
            a = a/Nbins
            expect_with.append(a)
            a = 0

    variance = Variance(expect)
    variance_with = Variance(expect_with)

    sigma.append(np.sqrt(variance))
    sigma_with.append(np.sqrt(variance_with))

fig = plt.figure()

ax = fig.add_subplot(111,xlabel='Nbins',ylabel='Standard deviation \sigma',xscale='log')

ax.plot(x,sigma,label=r'With $S_E$')
ax.plot(x,sigma_with,label=r'With $S_{eff}$')

ax.legend(loc='best')

fig.tight_layout()

fig.savefig('./plot/autocorrelation.pdf', bbox_inches='tight',format='pdf')
