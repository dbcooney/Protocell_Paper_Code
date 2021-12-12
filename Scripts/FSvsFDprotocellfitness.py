"""
Script used to generate Figure 4.4, illustrating the average protocell-level replication
rate / fitness of the population at steady state on both the fast-dimer (red) and
fast-slow (blue) edges of the simplex plotted as a function of the relative strength
of protocell-level competition $\lambda$. This figure highlights both a case in which
fast-slow competition can produce a better collective outcome than fast-dimer competition
for some intermediate values of $\lambda$ (left panel), as well as a case in which the
fast-dimer competition always produces an equal or better collective outcome
compared to that produced by fast-slow competition (right panel) for any $\lambda$.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import os

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


s = 1.
bF = 1. + s
bS = 1.
bD = (1. + s) / (2. + s)

theta = 2.

lamb_max = 300.
lamb_step = 0.1
lamb_range = np.arange(0.,lamb_max + lamb_step, lamb_step)



"""
Defining protocell-level replication rates and average collective replication rate
achieved at steady state for both the fast-slow and fast-dimer edges of the simplex.
"""

def GFS(x,eta):
	return x - eta * (x ** 2.)
	
def GFD(z,eta):
	return 0.5 * z * (1. - 0.5 * eta * z)
	
def GFSavg(lamb,eta, bF, bS, theta):
	if lamb * (1. - eta) < (bF - bS) * theta:
		return GFS(0.,eta)
	else:
		return  GFS(1.,eta) - ((bF - bS) * theta) / lamb
		
def GFDavg(lamb,eta, bF, bD, theta):
	if lamb * GFD(1.,eta) < (bF - bD) * theta:
		return GFD(0.,eta)
	else:
		return  GFD(1.,eta) - ((bF - bD) * theta) / lamb
		
GFSavgvec = np.vectorize(GFSavg)
GFDavgvec = np.vectorize(GFDavg)

plt.figure(1)

"""
Plotting collective outcomes on both edges of the simplex for a case in which fast-dimer
competition produces equal or better collective outcomes than the outcomes produced 
by fast-slow competition (when $\eta = 0.9$).
"""

eta = 0.9
		
plt.plot(lamb_range,GFSavgvec(lamb_range,eta,bF,bS,theta), lw = 6., color = 'b', label = r"$\langle G_{FS}(\cdot) \rangle_{f^{\lambda}_{\theta}}$")
plt.plot(lamb_range,GFDavgvec(lamb_range,eta,bF,bD,theta), lw = 6., color = 'r', label = r"$\langle G_{FD}(\cdot) \rangle_{g^{\lambda}_{\theta}}$")

plt.axhline(y = GFS(1.,eta), lw = 6., color = 'b', alpha = 0.8, ls = '--', label = r"$G_{FS}(1)$")
plt.axhline(y = GFD(1.,eta), lw = 6., color = 'r', alpha = 0.8, ls = '--', label = r"$G_{FD}(1)$")

plt.axis([0.,200.,0.,0.45])

plt.xlabel(r"Strength of Between-Protocell Selection $\lambda$", fontsize = 24., labelpad = 10.)
plt.ylabel(r"Average Protocell Fitness", fontsize = 24.)

plt.legend(loc = "upper center", ncol = 2, frameon = False, prop = {"size" : 22})

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.title(r"$\eta = 0.9$", fontsize = 24.)

plt.tight_layout()


script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)
plt.savefig(protocell_folder + "/Figures/FSvsFDprotocellfitnesseta0p9.png")

plt.figure(2)

"""
Plotting collective outcomes on both edges of the simplex for a case in fast-slow 
competition produces better collective outcomes than fast-diemr competition for
intermediate values of $\lambda$, but fast-dimer competition is more successful at
large values of $\lambda$ (when $\eta = 0.705$).
"""

eta = 0.705
	
plt.plot(lamb_range,GFDavgvec(lamb_range,eta,bF,bD,theta), lw = 6., color = 'r', label = r"$\langle G_{FD}(\cdot) \rangle_{g^{\lambda}_{\theta}}$")	
plt.plot(lamb_range,GFSavgvec(lamb_range,eta,bF,bS,theta), lw = 6., color = 'b', label = r"$\langle G_{FS}(\cdot) \rangle_{f^{\lambda}_{\theta}}$")

plt.axhline(y = GFD(1.,eta), lw = 6., color = 'r', alpha = 0.8, ls = '--', label = r"$G_{FD}(1)$")
plt.axhline(y = GFS(1.,eta), lw = 6., color = 'b', alpha = 0.8, ls = '--', label = r"$G_{FS}(1)$")


plt.axis([0.,100.,0.,0.5])

plt.xlabel(r"Strength of Between-Protocell Selection $\lambda$", fontsize = 24., labelpad = 10.)
plt.ylabel(r"Average Protocell Fitness", fontsize = 24.)

plt.legend(loc = "upper center", ncol = 2, frameon = False, prop = {"size" : 22})

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.title(r"$\eta = 0.705$", fontsize = 24.)

plt.tight_layout()


plt.savefig(protocell_folder + "/Figures/FSvsFDprotocellfitnesseta0p705.png")

plt.show()
		
	
