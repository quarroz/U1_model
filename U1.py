import numpy as np
import scipy.special as sp
import scipy.interpolate as si
from math import isfinite
import pickle
import sys

# This code follows the paper arXiv:1308.0233v1
# We apply the Lefschetz thimble method to one plaquette model.
# This code uses order 4 Runge-Kutta with adaptative stepsize
#
# The code demands an integer between 0 and 3 as an input
# to select the Nt we want to deal with in Nt_list

# We import the parameter written in param.py
import param

Nt_list = param.Nt_list
dr0 = param.dr0
dr_max = param.dr_max
N = param.N
precision = param.precision
f = param.f
beta_min = param.beta_min
beta_max = param.beta_max
beta_nb  = param.beta_nb
with_jacobian = param.with_jacobian
jac_evolve = param.jac_evolve

beta_list = np.linspace(beta_min,beta_max,beta_nb)

# Read which Nt we consider
j = int(sys.argv[1])
Nt = Nt_list[j]
epsilon = 1.0 /Nt
# Read which beta we consider
k = int(sys.argv[2])
runbeta = beta_list[k]

######################### Functions #####################################

def S(phi,beta=runbeta):
    """ Return S = -1j*beta*np.cos(phi) the action for one plaquette
    for given phi"""
    return -1j*beta*np.cos(phi)

def SG(eta,beta=runbeta):
    """ Return the gaussian action SG for a given complex phi"""
    return beta*eta**2/2.0

def S_eff(phi,jac,eta,beta=runbeta):
    return np.real(S(phi)) - np.log(abs(jac))

def step_RK4(phi,r,dr,beta=runbeta):
    """ Return the result of one step RK4 starting from pi,qi
    at position r with step dr """
    phi1 = 1.0/(r       )*np.conj(1j*beta*np.sin(phi         ))*dr
    phi2 = 1.0/(r+0.5*dr)*np.conj(1j*beta*np.sin(phi+0.5*phi1))*dr
    phi3 = 1.0/(r+0.5*dr)*np.conj(1j*beta*np.sin(phi+0.5*phi2))*dr
    phi4 = 1.0/(r+    dr)*np.conj(1j*beta*np.sin(phi+    phi3))*dr
    phi_next = phi + (phi1 + 2*phi2 + 2*phi3 + phi4)/6.0
    return phi_next

def RK4(phi,r,beta=runbeta):
    """ Runge-Kutta order 4 method with adaptative stepsize
    Return the last value p,q, real and imaginary part of the solution
    and the list of r and the corresponding hessian value."""
    dr = dr0
    r_list = [r]
    H = [hess_th(phi)]
    finite = True
    while r<1.0 and finite:
        # we make the calculation at dr and dr/2 and see if the diff d
        # between the two obtained value is below the precision
        phi_next = step_RK4(phi,r,dr,beta)
        phi_tmp = step_RK4(phi,r,dr*0.5,beta)
        phi_next_half = step_RK4(phi_tmp,r+dr*0.5,dr*0.5,beta)
        tmp = abs(phi_next - phi_next_half)
        if not(isfinite(np.imag(phi_next))) or np.real(phi_next)>np.pi/2:
            finite = False
            break
        if tmp == 0:
            dr = 10*dr
            continue
        if tmp<=precision:
            dr_old = dr
            dr = dr*(precision/tmp)**0.2  # Find the new stepsize
            if dr > dr_max:
                dr = dr_max
            r = r + dr_old
            if r+dr > 1:
                dr = 1-r
                phi = step_RK4(phi_next,r,dr)
                r = 1
                r_list.append(r)
                H.append(hess_th(phi))
            else:
                phi = phi_next
                r_list.append(r)
                H.append(hess_th(phi))
        else:
            dr = f*dr*(precision/tmp)**0.2
    return phi, r_list, H, finite

def hess_th(phi,beta=runbeta):
    """ Theoretical value of the Hessian evaluated in phi"""
    return 1j*beta*np.cos(phi)

def J_step(jac,H,r,dr,beta=runbeta):
    """ One step for the integration of the jacobian """
    j1 = 1.0/(r       )*np.conj(H(r       )*(jac        ))*dr
    j2 = 1.0/(r+dr/2.0)*np.conj(H(r+0.5*dr)*(jac+0.5*j1))*dr
    j3 = 1.0/(r+dr/2.0)*np.conj(H(r+0.5*dr)*(jac+0.5*j2))*dr
    j4 = 1.0/(r+dr)    *np.conj(H(r+    dr)*(jac+    j3))*dr
    jac_next = jac + 1.0/6.0*(j1+2.0*j2+2.0*j3+j4)
    return jac_next

def RK4_jac(j0,hess,r,beta=runbeta):
    """ Runge-Kutta order 4 method with adaptative stepsize for jacobian
    Return the last value p,q, real and imaginary part of the solution"""
    dr = dr0
    jac = j0
    while r<1.0:
        # we make the calculation at dr and dr/2 and see if the diff d
        # between the two obtained value is below the precision
        j_next = J_step(jac,hess,r,dr,runbeta)
        j_tmp = J_step(jac,hess,r,dr*0.5,beta)
        j_next_half = J_step(j_tmp,hess,r+dr*0.5,dr*0.5,beta)
        tmp = abs(j_next - j_next_half)
        if tmp == 0:
            dr = 10*dr
            continue
        if tmp<=precision:
            dr_old = dr
            dr = dr*(precision/tmp)**0.2  # Find the new stepsize
            r = r + dr_old
            if r+2*dr > 1:# Note the factor 2 otherwise H(r>1) is eval.
                dr = 1-r
                jac = J_step(j_next,hess,r,dr,beta)
                r = 1
            else:
                jac = j_next
        else:
            dr = f*dr*(precision/tmp)**0.2
    jac = jac * epsilon**beta # To obtain J^phi_eta
    return jac

def jac_th(phi,eta,beta=runbeta):
    """ Return the theoretical value of the jacobian """
    if Nt==1:
        return 1+1j*0
    else:
        return -1j*np.conj(np.sin(phi))/eta

def residual(field,eta,beta=runbeta):
    """ Take as input a complex field list and an eta list and return
    the value of the residual which is a complex number """
    result = jac_th(field,eta,beta)*np.exp(-S(field))
    return result

def Metropolis(phi_old, eta_old, jac_old,
               phi_new, eta_new, jac_new, beta=runbeta):
    if with_jacobian:
        value = np.exp(-S_eff(phi_new, jac_new, eta_new)\
                   +S_eff(phi_old, jac_old, eta_old)\
                   +SG(eta_new,beta)-SG(eta_old,beta))
    else:
        value = np.exp(-np.real(S(phi_new))+np.real(S(phi_old))\
                 +SG(eta_new,beta)-SG(eta_old,beta))
    Proba_accept = min(1,value)
    x = np.random.random() # Sample a random number between [0,1)
    return x<= Proba_accept

def Fa(phi,jac,beta=runbeta):
    if not(with_jacobian):
        return np.exp(-1j*np.imag(S(phi,beta)))*epsilon**beta*jac
    else:
        return np.exp(-1j*np.imag(S(phi,beta)))*epsilon**beta\
                *np.exp(1j*np.angle(jac))

def expect_value(phi,jac):
    num = 1j*np.imag(np.exp(1j*phi)*Fa(phi,jac))
    denom = np.real(Fa(phi,jac))
    return num, denom

######################## Initialize observable ##########################

accept_p = []
accept_q = []
x4 = []
y4 = []
jac_th_list = []
jac_list = []
numerator = 0
Z = 0
convergence = []

########################### Parameters ##################################

w = 1/np.sqrt(2)*(1-1j)
phi_c = 0 +1j*0
j0 = w

######################### Main program ##################################

phi_old = phi_c # The first field is initialize for Metropolis
eta_old = 0+0*1j
jac_old = w
i = 0
j = 0

print('Starting process for N={0},beta={1} and Nt={2}'.format(N,runbeta,Nt))

while i<N: #Loop over iteration
    phi = 4+4*1j # Just there to enter in the while loop
    finite = False
    while abs(np.real(phi))>np.pi/2 or not(finite):
        # Check that the initial field is finite
        eta = np.random.normal(0,1/np.sqrt(runbeta))
        phi = phi_c + w*epsilon**runbeta*eta
        if np.real(phi)>np.pi/2: # It is not necessary to go further ...
            continue
        phi, r_list, H, finite = RK4(phi,epsilon)
    # Jacobian evolution
    if jac_evolve:
        if Nt!=1: # Do evolution of jacobian with interpolation
            hess = si.interp1d(r_list,H,'slinear',assume_sorted=True,
                           bounds_error=False,fill_value='extrapolate')
            jac = RK4_jac(j0,hess,epsilon)
        else: # Gaussian case
            jac = j0
    else:
        if Nt!=1:
            jac = jac_th(phi,eta,runbeta)
        else:
            jac = j0

    i+=1
    # Metropolis test
    if Metropolis(phi_old, eta_old, jac_old, phi, eta, jac):
        phi_old = phi
        eta_old = eta
    else:
        j+=1
        continue
    # Print
    #print('phi',phi)
    #print('eta',eta)
    #print('jac',jac)
    #if Nt!=1:
    #    print('jac_th',jac_th(phi,eta,runbeta))
    #else:
    #    print('jac_th',j0)

    # Fill observable
    accept_p.append(np.real(phi))
    accept_q.append(np.imag(phi))
    res = residual(phi,eta)
    x4.append(abs(res))
    y4.append(np.cos(np.angle(res)))
    jac_th_list.append(jac_th(phi,eta,runbeta))
    jac_list.append(jac)

    # Expectation value calculation
    num, denom = expect_value(phi,jac)
    numerator = (numerator*i+num)/(i+1)
    Z = (Z*i+denom)/(i+1)
    convergence.append(np.imag(numerator/Z))

print(j)
expect = convergence[-1]
print('Nt={0}, beta={1}, result={2}'.format(Nt,runbeta,expect))

########################### Write data in a file ########################
with open("data/accept_p_Nt_{0}_B_{1}.txt".format(Nt,k),"wb") as fp:
    pickle.dump(accept_p,fp)
with open("data/accept_q_Nt_{0}_B_{1}.txt".format(Nt,k),"wb") as fp:
    pickle.dump(accept_q,fp)
with open("data/expect_Nt_{0}_B_{1}.txt".format(Nt,k),"wb") as fp:
    pickle.dump(expect,fp)
with open("data/x4_Nt_{0}_B_{1}.txt".format(Nt,k),"wb") as fp:
    pickle.dump(x4,fp)
with open("data/y4_Nt_{0}_B_{1}.txt".format(Nt,k),"wb") as fp:
    pickle.dump(y4,fp)
with open("data/jac_th_Nt_{0}_B_{1}.txt".format(Nt,k),"wb") as fp:
    pickle.dump(jac_th_list,fp)
with open("data/jac_Nt_{0}_B_{1}.txt".format(Nt,k),"wb") as fp:
    pickle.dump(jac_list,fp)
with open("data/convergence_Nt_{0}_B_{1}_jac_evolve_{2}.txt".format(Nt,k,with_jacobian),"wb") as fp:
    pickle.dump(convergence,fp)
