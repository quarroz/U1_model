import numpy as np

################################ Parameters ############################

# Parameters
Nt_list=[1,2,20,200]
default_beta = 3
dr0 = 1e-10
N_step = 1e6
dr_max = 1# 0.01
N = 8000
precision = 1e-8     # precision for adaptative stepsize
f = 0.95             # avoid infinite loop when reducing stepsize
choosen_beta = 8
T_list = [0.0,0.1,0.2,0.3]

real = np.linspace(-np.pi/2.0,np.pi/2.0)
imag = np.linspace(-np.pi,np.pi)
beta_min = 0.4
beta_max = 2.4
beta_nb  = 11

marker_list = ["v",".","^","s"]
color_list = ["red","darkorange","forestgreen","royalblue"]

with_jacobian = False
jac_evolve = False
