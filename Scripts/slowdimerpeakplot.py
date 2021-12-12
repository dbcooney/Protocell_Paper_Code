"""
Script used to generate Figure 4.8, which provides a comparison between the protocell
composition with maximal protocell-level replication rate and the modal composition 
achieved at steady state for dimorphic competition on the slow-dimer edge of the 
simplex in the limit of strong between-group competition (as $\lambda \to \infty$).
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spi
import os

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


"""
Formulas characterizing composition with maximal protocell-level replication rate and
modal composition of the population at steady states in the limit of strong 
between-protocell competition.
"""

		
def peak_group(eta):
	if eta <= 2.0/3.0:
		return 0.
	else:
		return 3.0 - (2.0 / eta)
		
def max_payoff_group(eta):
	if eta <= 0.5:
		return 0.
	else:
		return 2.0 - (1.0 / eta)
		


max_dimer_vec = np.vectorize(max_payoff_group)		
peak_dimer_vec = np.vectorize(peak_group)	

	
eta_step = 0.01
eta_range = np.arange(0.0,1.0+eta_step,eta_step)
	


plt.axvline(x = 0.5, lw = 5., ls = '--', alpha = 0.8, color = 'gray')
plt.axvline(x = 2.0/3.0, lw = 5., ls = '--', alpha = 0.8, color = 'gray')

plt.plot(eta_range,max_dimer_vec(eta_range), lw = 6., color = 'b', label = r"$z^*_{SD}$: Optimal Protocell Composition ")
plt.plot(eta_range,peak_dimer_vec(eta_range), lw = 6., color = 'g', ls = '-', label = r"$\hat{z}_{\infty}$: Modal Composition")


plt.axis([0.0,1.0,-0.01,1.0])	

plt.xlabel(r"Complementarity Parameter $\eta$", fontsize = 24., labelpad = 10.)
plt.ylabel(r"Fraction of Dimers $z$", fontsize = 24.)

plt.legend(loc = "upper left", prop = {"size":18})

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.tight_layout()



script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)
plt.savefig(protocell_folder + '/Figures/slowdimerpeakplot.png')


plt.show()