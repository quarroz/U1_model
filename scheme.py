import matplotlib.pyplot as plt
import numpy as np

#LaTeX
plt.rc('text',usetex=True)
plt.rc('font',family='serif')

scheme = plt.figure()

def exact_thimble(x):
    result = []
    for i in x:
        if i>=0:
            result.append(-np.arccosh(1.0/np.cos(i)))
        else:
            result.append(np.arccosh(1.0/np.cos(i)))
    return result

def exact_thimble2(x):
    result = []
    for i in x:
        if i < np.pi:
            result.append(-np.arccosh(-1.0/np.cos(i)))
        else:
            result.append(np.arccosh(-1.0/np.cos(i)))
    return result

def multiple_formatter(denominator=2, number=np.pi, latex='\pi'):
    def gcd(a, b):
        while b:
            a, b = b, a%b
        return a
    def _multiple_formatter(x, pos):
        den = denominator
        num = np.int(np.rint(den*x/number))
        com = gcd(num,den)
        (num,den) = (int(num/com),int(den/com))
        if den==1:
            if num==0:
                return r'$0$'
            if num==1:
                return r'$%s$'%latex
            elif num==-1:
                return r'$-%s$'%latex
            else:
                return r'$%s%s$'%(num,latex)
        else:
            if num==1:
                return r'$\frac{%s}{%s}$'%(latex,den)
            elif num==-1:
                return r'$\frac{-%s}{%s}$'%(latex,den)
            else:
                return r'$\frac{%s%s}{%s}$'%(num,latex,den)
    return _multiple_formatter

class Multiple:
    def __init__(self, denominator=2, number=np.pi, latex='\pi'):
        self.denominator = denominator
        self.number = number
        self.latex = latex

    def locator(self):
        return plt.MultipleLocator(self.number / self.denominator)

    def formatter(self):
        return plt.FuncFormatter(multiple_formatter(self.denominator, self.number, self.latex))

ax = scheme.add_subplot(111)
ax.set(xlim=[-2,5])
ax.set(ylim=[-np.pi,np.pi])
ax.set_aspect(1.0)
ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 2))
ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))

plt.axvline(x=np.pi/2,linestyle='dashed',color='red')
plt.axvline(x=-np.pi/2,linestyle='dashed',color='red')
plt.axvline(x=3*np.pi/2,linestyle='dashed',color='red')
plt.axvline(x=-np.pi/2,ymin=0.5,color='green')
plt.axvline(x=3/2*np.pi,ymin=0.5,color='green')
plt.axhline(y=0,xmin=(-np.pi/2+2)/7,xmax=(3/2*np.pi+2)/7,color='blue')

x=np.linspace(-np.pi/2,np.pi/2,10000)
ax.plot(x,exact_thimble(x),color='green')

y=np.linspace(np.pi/2,3*np.pi/2,10000)
ax.plot(y,exact_thimble2(y),color='green')

s = 'scheme_deformation'

scheme.savefig(s+'.pdf', bbox_inches='tight',format='pdf')
