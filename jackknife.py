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
        result+=el
    return result/N

def Jackknife(ensemble):
    final = []
    N = len(ensemble)
    for i in range(len(ensemble)):
        result = 0
        for j in range(len(ensemble)):
            if j!=i:
                result += ensemble[j]
        result = result/(N-1)
        final.append(result)
    return final

def Variance(ensemble):
    result = 0
    N = len(ensemble)
    moy = Moyenne(ensemble)
    for el in ensemble:
        result += (el - moy)**2
    return result*(N-1)/N


######################### Reading the data #############################
#for i in range(1000):
#    with open("data/expectation_{0}_False.txt".format(i),"rb") as fp:
#        expect_full.append(pickle.load(fp))
#    with open("data/expectation_{0}_True.txt".format(i),"rb") as fp:
#        expect_with_full.append(pickle.load(fp))

expect_full = [1,2,3,4,5,6,7,8,9,10] 

beta = 1.0
Nt = 20
N=len(expect_full)
text = r'Method: Jackknife\\'+r'Nt={0} - $\beta$={1}\\'.format(Nt,beta) + r'N={0}\\'.format(N)
props = dict(boxstyle='round',facecolor='whitesmoke',alpha=0.7)

expect = Jackknife(expect_full)
expect_with = Jackknife(expect_full)

print('Initial length: ', len(expect_full))
print('Final length: ', len(expect))

mean = Moyenne(expect)
mean_with = Moyenne(expect_with)

variance = Variance(expect)
variance_with = Variance(expect_with)

sigma = np.sqrt(variance)
sigma_with = np.sqrt(variance_with)

print(sigma)

result = r'<$\langle exp(i\phi)\rangle = {:.5f} \pm {:.5f}$'.format(mean,sigma)
result_with = r'$\langle exp(i\phi) \rangle = {:.5f} \pm {:.5f}$'.format(mean_with,sigma_with)

print(r'Avec $S_E$ :'+result)
print(r'Avec $S_{eff}$ :'+result_with)

error = plt.figure()
error_with = plt.figure()

x = np.linspace(min(expect),max(expect),100)
x_with = np.linspace(min(expect_with),max(expect_with),100)

ax = error.add_subplot(111,xlabel=r'Expectation value of $e^{i\phi}$',
        ylabel='Occurrence')
#ax.set(xlim=[min(expect),max(expect)])
ax_with = error_with.add_subplot(111,xlabel=r'Expectation value of $e^{i\phi}$',
        ylabel='Occurrence')
#error_with.set(xlim=[min(expect_with),max(expect_with)])

print('Here', len(expect))

ax.hist(expect,label=r'Using $S_E$')
ax_with.hist(expect_with, label=r'Using $S_{eff}$')

#ax.plot(x,mlab.normpdf(x,mean,sigma))
#ax_with.plot(x_with,mlab.normpdf(x_with,mean_with,sigma_with))

error.text(0.6,0.4,text+result,transform=ax.transAxes,fontsize=14, verticalalignment='bottom',bbox=props)

error_with.text(0.6,0.4,text+result_with,transform=ax_with.transAxes,fontsize=14, verticalalignment='bottom',bbox=props)

ax.legend(loc='best')
ax_with.legend(loc='best')

error.tight_layout()
error_with.tight_layout()

error.savefig('./plot/jackknife.pdf', bbox_inches='tight',format='pdf')
error_with.savefig('./plot/jackknife_with.pdf', bbox_inches='tight',format='pdf')

with open("data.txt","w") as f:
    for el in expect_full:
        f.write(str(el)+"\n")
