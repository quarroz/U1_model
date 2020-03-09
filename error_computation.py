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
    return result * 1/(N-1)


######################### Reading the data #############################
for i in range(1000):
    with open("data/expectation_{0}_False.txt".format(i),"rb") as fp:
        expect_full.append(pickle.load(fp))
    with open("data/expectation_{0}_True.txt".format(i),"rb") as fp:
        expect_with_full.append(pickle.load(fp))

beta = 1.0
Nt = 20
Nbins=int(sys.argv[1])
N=len(expect_full)
text = r'Method: binning\\'+r'Nt={0} - $\beta$={1}\\'.format(Nt,beta) + r'Nbins={0} - N={1}\\'.format(Nbins,N)
props = dict(boxstyle='round',facecolor='whitesmoke',alpha=0.7)

a=0
for i in range(1,len(expect_full)+1):
    if i%Nbins!=0:
        a+=expect_full[i-1]
    else:
        a+=expect_full[i-1]
        a = a/Nbins
        expect.append(a)
        a = 0

a=0
for i in range(1,len(expect_with_full)+1):
    if i%Nbins!=0:
        a+=expect_with_full[i-1]
    else:
        a+=expect_full[i-1]
        a = a/Nbins
        expect_with.append(a)
        a = 0

mean = Moyenne(expect)
mean_with = Moyenne(expect_with)

print(mean_with)
print(np.mean(expect_with))

variance = Variance(expect)
variance_with = Variance(expect_with)

sigma = np.sqrt(variance)
sigma_with = np.sqrt(variance_with)

result = r'<$\langle exp(i\phi)\rangle = {:.3f} \pm {:.3f}$'.format(mean,sigma)
result_with = r'$\langle exp(i\phi) \rangle = {:.3f} \pm {:.3f}$'.format(mean_with,sigma_with)

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

ax.hist(expect,label=r'Using $S_E$',density=True)
ax_with.hist(expect_with, label=r'Using $S_{eff}$',density=True)

ax.plot(x,mlab.normpdf(x,mean,sigma))
ax_with.plot(x_with,mlab.normpdf(x_with,mean_with,sigma_with))

error.text(0.6,0.4,text+result,transform=ax.transAxes,fontsize=14, verticalalignment='bottom',bbox=props)

error_with.text(0.6,0.4,text+result_with,transform=ax_with.transAxes,fontsize=14, verticalalignment='bottom',bbox=props)

ax.legend(loc='best')
ax_with.legend(loc='best')

error.tight_layout()
error_with.tight_layout()

error.savefig('./plot/error_histo.pdf', bbox_inches='tight',format='pdf')
error_with.savefig('./plot/error_with_histo.pdf', bbox_inches='tight',format='pdf')
