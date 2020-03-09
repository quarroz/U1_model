import matplotlib.pyplot as plt
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
accept_p = []
accept_q = []
expect = []
field_r = []
field_i = []
x4 = []
y4 = []
convergence = []
convergence_with = []

######################### Reading the data #############################
for Nt in Nt_list:
    with open("data/accept_p_Nt_{0}_B_{1}.txt".format(Nt,default_beta),"rb") as fp:
        accept_p.append(pickle.load(fp))
    with open("data/accept_q_Nt_{0}_B_{1}.txt".format(Nt,default_beta),"rb") as fp:
        accept_q.append(pickle.load(fp))
    with open("data/accept_p_Nt_{0}_B_{1}.txt".format(Nt,choosen_beta),\
              "rb") as fp:
        field_r.append(pickle.load(fp))
    with open("data/accept_q_Nt_{0}_B_{1}.txt".format(Nt,choosen_beta),\
              "rb") as fp:
        field_i.append(pickle.load(fp))
    with open("data/x4_Nt_{0}_B_{1}.txt".format(Nt,default_beta),"rb") as fp:
        x4.append(pickle.load(fp))
    with open("data/y4_Nt_{0}_B_{1}.txt".format(Nt,default_beta),"rb") as fp:
        y4.append(pickle.load(fp))
#    with open("data/jac_th_Nt_{0}_B_{1}.txt".format(Nt,default_beta),"rb") as fp:
#        jac_th_list.append(pickle.load(fp))
#    with open("data/jac_Nt_{0}_B_{1}.txt".format(Nt,default_beta),"rb") as fp:
#        jac.append(pickle.load(fp))
    with open("data/convergence_Nt_{0}_B_{1}_jac_evolve_False.txt".format(Nt,default_beta),"rb") as fp:
        convergence.append(pickle.load(fp))
    with open("data/convergence_Nt_{0}_B_{1}_jac_evolve_True.txt".format(Nt,default_beta),"rb") as fp:
        convergence_with.append(pickle.load(fp))
    list_Nt = []
    for i in range(beta_nb):
        with open("data/expect_Nt_{0}_B_{1}.txt".format(Nt,i),\
                  "rb") as fp:
            list_Nt.append(pickle.load(fp))
    expect.append(list_Nt)

default_beta = beta[default_beta]
choosen_beta = beta[choosen_beta]
########################## Plot function ###############################

def exact_thimble(x):
    result = []
    for i in x:
        if i>=0:
            result.append(-np.arccosh(1.0/np.cos(i)))
        else:
            result.append(np.arccosh(1.0/np.cos(i)))
    return result

def expect_th(beta=default_beta):
    # Return the imaginary part of the theoretical expect. value
    result = []
    for i in beta:
        result.append(sp.j1(i)/sp.j0(i))
    return result


########################## Plot section ################################
thimble_default = plt.figure()
thimble_choosen = plt.figure()
rp = plt.figure()
expect_fig = plt.figure()
#convergence_fig1 = plt.figure()
#convergence_fig2 = plt.figure()

# Setup + theoretical values

# First plot - Two distributions
ax = thimble_default.add_subplot(111,xlabel='Real part of the field',\
                ylabel='Imaginary part of the field',\
                title=r'Sampled field at $\beta={:04.2f}$'.format(default_beta))
ax.set(xlim=[-np.pi/2.0,np.pi/2.0],ylim=[-np.pi,np.pi])
# Plot the theorical value
ax.plot(real,exact_thimble(real),color='black',alpha=0.2,\
        label='exact',linewidth='2')

ax2 = thimble_choosen.add_subplot(111,xlabel='Real part of the field',\
                ylabel='Imaginary part of the field',\
                title=r'Sampled field at $\beta={:04.2f}$'.format(choosen_beta))
ax2.set(xlim=[-np.pi/2.0,np.pi/2.0],ylim=[-np.pi,np.pi])
# Plot the theorical value
ax2.plot(real,exact_thimble(real),color='black',alpha=0.2,\
        label='exact',linewidth='2')

# Expectation value
title = r'Expect. value of $e^{i\phi}$ wrt $\beta$'
ax3 = expect_fig.add_subplot(111,title=title,ylabel=r'$\Im e^{i\phi}$',\
                  xlabel=r'$\beta$',yscale='log')
ax3.set(xlim=[0,2.5]) # 2.4 value of first non-trivial zero of J0
# Plot theoretical value
beta_th = np.linspace(0.1,2.4,500)
ax3.plot(beta_th,expect_th(beta_th),color='black',alpha=0.2,\
         linewidth='2.5', linestyle='-',label='exact')

# Residual phase
ax4 = rp.add_subplot(111,xlabel=r'$|J^\phi_\eta e^{-S}|$',\
        ylabel=r'$cos(arg\{J^\phi_\eta e^{-S}\})$',xscale='log')

# Convergence
#N_list = np.linspace(1,N,N)
#ax5 = convergence_fig1.add_subplot(111, xlabel=r'Number of element in statistical ensemble', ylabel=r'Relative error on $\Im<e^{i\phi}>$ $[\%]$')
#ax5.set(ylim=[0,0.2])
#ax6 = convergence_fig2.add_subplot(111, xlabel=r'Number of element in statistical ensemble', ylabel=r'Relative error on $<e^{i\phi}>$ $[\%]$')
#ax6.set(ylim=[0,0.2])
#######

# Data treatment

#convergence = (convergence - sp.j1(default_beta)/sp.j0(default_beta))/(sp.j1(default_beta)/sp.j0(default_beta))*100

# Plot the data
for i in range(len(Nt_list)):
    Nt = Nt_list[i]
    # Esthetic
    color = color_list[i]
    marker = marker_list[i]
    if Nt!=1:
        label=r"$N_\tau={0}$".format(Nt)
    else:
        label='gaussian'
    alpha = 0.5

    # Plot the field configuration at default beta
    ax.scatter(accept_p[i],accept_q[i],marker=marker,label=label,\
               edgecolors=color,facecolors=color,s=22, linewidths=0.3,alpha=alpha-0.2)

    # Plot expectation value of exp(i*phi)
    ax3.scatter(beta,expect[i],marker=marker,edgecolors=color,\
                    label=label,s=28.5,facecolors=color,linewidths=0.3,alpha=alpha+0.3)

    # Plot the field configuration at choosen beta
    ax2.scatter(field_r[i],field_i[i], marker=marker, label=label,\
            edgecolors=color,facecolors=color,s=22, linewidths=0.3,alpha=alpha-0.2)

    # Plot the residual
    if i==0:
        pass
    else:
        ax4.scatter(x4[i],y4[i],marker=marker,label=label, edgecolors=color,\
               facecolor=color,s=15.5, linewidth=0.7,alpha=1)

#    if i==1:
#        # Plot the convergence with jacobian at default beta
#        ax5.scatter(N_list,convergence_with[i],marker='.',
#            label=label+r' with $S_{eff}$',edgecolors=color_list[1],facecolors=color_list[1],s=18.5, linewidths=0.3,alpha=alpha)
#
#        # Plot plot convergence without jacobian at default beta
#        ax5.scatter(N_list,convergence[i],marker='.',
#            label=label,edgecolors=color_list[0],facecolors=color_list[0],
#                s=18.5, linewidths=0.3,alpha=alpha)
#
#    if i==3:
#        # Plot the convergence with jacobian at default beta
#        ax6.scatter(N_list,convergence_with[i],marker='.',
#            label=label+r' with $S_{eff}$',edgecolors=color_list[1],facecolors=color_list[1], s=18.5, linewidths=0.3, alpha=alpha)
#        # Plot plot convergence without jacobian at default beta
#        ax6.scatter(N_list,convergence[i],marker='.',
#            label=label,edgecolors=color_list[0],facecolors=color_list[0],
#                s=18.5, linewidths=0.3, alpha=alpha)

# Zoom on data


# Place the legend 
ax.legend(loc='best',prop={'size':16})
ax2.legend(loc='best',prop={'size':16})
ax3.legend(loc='best',prop={'size':18})
ax4.legend(loc='best',prop={'size':17})
#ax5.legend(loc='best',prop={'size':12})
#ax6.legend(loc='best',prop={'size':12})

# Name of the figure
s1 = 'plot/Nt={0}_{1}_{2}_{3}_N={4}_field'.format(Nt_list[0],\
                                    Nt_list[1],Nt_list[2],Nt_list[3],N)
s12 = 'plot/Nt={0}_{1}_{2}_{3}_N={4}_field2'.format(Nt_list[0],\
                                    Nt_list[1],Nt_list[2],Nt_list[3],N)
s3 = 'plot/Nt={0}_{1}_{2}_{3}_N={4}_expect'.format(Nt_list[0],\
                                    Nt_list[1],Nt_list[2],Nt_list[3],N)
s2 = 'plot/Nt={0}_{1}_{2}_{3}_N={4}_rp'.format(Nt_list[0],\
                                    Nt_list[1],Nt_list[2],Nt_list[3],N)
s4 = 'plot/Nt={0}_{1}_{2}_{3}_N={4}_convergence'.format(Nt_list[0],\
                                    Nt_list[1],Nt_list[2],Nt_list[3],N)
s5 = 'plot/Nt={0}_{1}_{2}_{3}_N={4}_convergence2'.format(Nt_list[0],\
                                    Nt_list[1],Nt_list[2],Nt_list[3],N)

thimble_default.tight_layout()
thimble_choosen.tight_layout()
rp.tight_layout()
expect_fig.tight_layout()
#convergence_fig1.tight_layout()
#convergence_fig2.tight_layout()

# Choose a type for saving the figure

#plt.savefig(s1+'.png', bbox_inches='tight',format='png')
thimble_default.savefig(s1+'.pdf', bbox_inches='tight',format='pdf')
#plt.savefig(s1+'.svg', bbox_inches='tight',format='svg')

#plt.savefig(s12+'.png', bbox_inches='tight',format='png')
thimble_choosen.savefig(s12+'.pdf', bbox_inches='tight',format='pdf')
#plt.savefig(s12+'.svg', bbox_inches='tight',format='svg')

#plt.savefig(s2+'.png', bbox_inches='tight',format='png')
rp.savefig(s2+'.pdf', bbox_inches='tight',format='pdf')
#plt.savefig(s2+'.svg', bbox_inches='tight',format='svg')

#plt.savefig(s3+'.png', bbox_inches='tight',format='png')
expect_fig.savefig(s3+'.pdf', bbox_inches='tight',format='pdf')
#plt.savefig(s3+'.svg', bbox_inches='tight',format='svg')

#plt.savefig(s4+'.png', bbox_inches='tight',format='png')
#convergence_fig1.savefig(s4+'.pdf', bbox_inches='tight',format='pdf')
#plt.savefig(s4+'.svg', bbox_inches='tight',format='svg')

#plt.savefig(s5+'.png', bbox_inches='tight',format='png')
#convergence_fig2.savefig(s5+'.pdf', bbox_inches='tight',format='pdf')
#plt.savefig(s5+'.svg', bbox_inches='tight',format='svg')
