"""
Script used to generate Figure 4.6, illustrating the parameter regimes in which fast-dimer 
or fast-slow dimorphic competition produces better outcomes via multilevel selection over
some range of the relative selection strength $\lambda$. The two parameters studied here
are the gene-level advantage for fast replicators $s$ and the complementarity parameter 
eta, and the regions in the figure correspond to regimes in which fast-dimer or fast-slow
competition features a lower threshold selection strength $\lambda^*$ to achieve
coexistence of the fast and slow genes and the regions in which fast-dimer or fast-slow
competition produces greater collective outcomes in the limit of infinite strength of
between-protocell competition.
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

s_max = 4.
eta_step = 0.001
eta_range = np.arange(0.0,1.0+eta_step,eta_step)
eta_23_range = np.arange(0.5,2.0/3.,eta_step)
eta_3p_range = np.arange(2.0/3.0 + eta_step, 1.0+ eta_step, eta_step)
print np.size(eta_range)


"""
Defining boundary between region in which fast-slow competition (orange) and in which
fast-dimer competition (red) has lower threshold selection strength.
""" 

def eta_thresh(eta):
	if eta <= 2.0/3.0:
		return -1.
	else:
		denom = (3.0* eta - 2.0)
		sqrt_term = ((3.0 * eta - 2.0)**2.0) - 4.0 * (eta - 1.0) * (3.0 * eta - 2.0)
		num = -(3.0 * eta - 2.0) + np.sqrt(sqrt_term)
		return num / denom
		

"""
Plotting the three regions of the relative collective benefits for fast-slow and 
fast-dimer competition. The regions depicted are those in which the slow-dimer
competition features both a lower threshold to achieve coexistence and greater collective 
outcomes in the limit of large between-protocell competition (yellow), in which the 
fast-slow competition has a lower threshold but fast-dimer competition produces greater
collective fitness for strong between-protocell competition (orange), and in which
fast-dimer competition featuress better outcomes across the range of relative selection
strengths (red).
"""

eta_thresh_vec = np.vectorize(eta_thresh)

thresh_holder = eta_thresh_vec(eta_range)
thresh_23_holder = eta_thresh_vec(eta_23_range)
thresh_3p_holder = eta_thresh_vec(eta_3p_range)
		
plt.plot(eta_range,eta_thresh_vec(eta_range), lw = 3., ls = '--', color = 'k')

plt.axvline(x = 2.0/3.0, lw = 3., ls = '--', color = 'k')


s_max_holder = [s_max] * np.size(eta_range)
s_max_23_holder = [s_max] * np.size(eta_23_range)
s_max_3p_holder = [s_max] * np.size(eta_3p_range)
zero_23_holder = [0.0] * np.size(eta_23_range)
zero_3p_holder = [0.0] * np.size(eta_3p_range)
plt.fill_between(eta_range,thresh_holder,s_max_holder, color = plt.cm.YlOrRd(0.1), alpha = 0.6)
plt.fill_between(eta_3p_range,thresh_3p_holder,zero_3p_holder , color = plt.cm.YlOrRd(0.5), alpha = 0.6)
plt.fill_between(eta_3p_range,s_max_3p_holder,thresh_3p_holder , color = plt.cm.YlOrRd(0.9), alpha = 0.6)

plt.annotate(r"$\lambda^*_{FD} > \lambda^*_{FS} \\ G_{FD}(1) < G_{FS}(1)$", xy = (0.15,2.9), fontsize = 18.)


plt.annotate(r"$\lambda^*_{FD} > \lambda^*_{FS} \\ G_{FD}(1) > G_{FS}(1)$",
            xy=(0.72, .53), xycoords='data',
            xytext=(0.375, 0.9), textcoords='data', fontsize = 18.,
            arrowprops=dict(arrowstyle="->", linewidth = 3))


plt.annotate(r"$\lambda^*_{FD} < \lambda^*_{FS} \\ G_{FD}(1) > G_{FS}(1)$ ", xy = (0.73,2.9), fontsize = 18.)





plt.axis([0.,1.,0.,s_max])


plt.xlabel(r"Complementarity Parameter $\eta$", fontsize = 24., labelpad = 10.)
plt.ylabel(r"Fast Replicator Advantage $s$", fontsize = 24.)

plt.tight_layout()


script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)
plt.savefig(protocell_folder + '/Figures/slowdimerregimes.png')

plt.show()
	