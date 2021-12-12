"""
Script used to generate Figure 4.5, which provides a comparison between the average
protocell-level replication rate (left panel) and the modal fraction of slow gene 
achieved at steady state (right panel) in the limit of infinite relative strength of
between-protocell competition ($\lambda \to \infty$), plotted as a function of the
complementarity parameter ($\eta$).
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


eta_step = 0.001
eta_range = np.arange(0.0,1.0+eta_step,eta_step)


"""
Defining collective protocell-level reproduction rate in the large-$\lambda$ limit
for the fast-dimer and fast-slow edges of the simplex.
"""

def GFD1(eta):
	return 0.5 * (1.0 - 0.5 * eta)
	
def GFS1(eta):
	return 1.0 - eta
	
"""
Defining modal fraction of slow replicators at steady state on the fast-slow edge of the
simplex in the limit as $\lambda \to \infty$. For the fast-dimer edge, the modal outcome
in this limit is the all-dimer composition, so the fraction of slow replicators in this 
case in 0.5 for any value of the complementarity parameter $\eta$. 
"""
	
def peakFS(eta):
	if eta < 0.5:
		return 1.0
	else:
		return (1.0 / eta) - 1.0
		
peakFS_vec = np.vectorize(peakFS)
	


plt.figure(1)

"""
Plotting comparison of collective payoffs as function of the complementarity parameter
$\eta$.
"""

plt.plot(eta_range,GFD1(eta_range), lw = 6., ls = '-', color = 'r', label = r"Fast-Dimer Edge")
plt.plot(eta_range,GFS1(eta_range), lw = 6., ls = '--', color = 'b', label = r"Fast-Slow Edge")

plt.axvline(x = 0.5, lw = 6., color = 'k', ls = '--', alpha = 0.8)
plt.axvline(x = 2.0/3.0, lw = 6., color = 'k', ls = '--', alpha = 0.8)


plt.axis([0.,1.0,0.0,1.01])


plt.ylabel(r"Maximal Steady State Payoff ($G_{\cdot}(1)$)", fontsize = 20.)
plt.xlabel(r"Complementarity Parameter $\eta$", fontsize = 25., labelpad = 20.)

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.legend(loc = "lower left")

plt.tight_layout()


script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)
plt.savefig(protocell_folder + '/Figures/fsfdpayoffcomparison.png')

plt.figure(2)

"""
Plotting comparison of modal fractions of slow genes as function of the complementarity
parameter $\eta$.
"""

plt.plot(eta_range,[0.5] * np.size(eta_range), lw = 6., ls = '-', color = 'r', label = r"Fast-Dimer Edge")
plt.plot(eta_range,peakFS_vec(eta_range), lw = 6., ls = '--', color = 'b', label = r"Fast-Slow Edge")



plt.axvline(x = 0.5, lw = 6., color = 'k', ls = '--', alpha = 0.8)
plt.axvline(x = 2.0/3.0, lw = 6., color = 'k', ls = '--', alpha = 0.8)

plt.xlabel(r"Complementarity Parameter $\eta$", fontsize = 25., labelpad = 20.)
plt.ylabel(r"Slow Types in Modal Group ($\lambda \to \infty$)", fontsize = 20.)

plt.axis([0.,1.0,0.0,1.01])

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.legend(loc = "lower left")

plt.tight_layout()

plt.savefig(protocell_folder + '/Figures/fsfdpeakcomparison.png')

plt.show()