"""
Script used to generate Figure 2.3, which provides a comparison between the protocell
composition with maximal protocell-level replication rate and the modal composition 
achieved at steady state for the baseline protocell model (on the fast-slow edge
of the simplex) in the limit of strong between-group competition 
(as $\lambda \to \infty$).
"""


import matplotlib.pyplot as plt
import numpy as np
import os



from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)




eta = 0.7
alpha = 2.
length = 1000
space_step = 0.025

lamb_list = []
peak1 = []
peak2 = []
peak_list = []
max_fit_list = int(length / space_step) * [1. / (2. * eta)]


lamb = 200

def steady_state(x):
	return (x ** (lamb * (1. - eta) - alpha - 1)) * ((1. - x)**(alpha - 1.)) * np.exp(- lamb * eta * x)
	
steady_state_vec = np.vectorize(steady_state)
current_steady_state = steady_state(np.arange(0.0,1.0,0.01))
current_steady_state = current_steady_state / np.sum(current_steady_state)
best_group = np.argmax(current_steady_state)
print best_group
	
	




eta_range = np.arange(0.,1.,0.01)

"""
Introducing formulas for maximum protocell-level replication rate and modal composition
at steady state as $\lambda \to \infty$, characterized as a function of the 
complementarity parameter $\eta$.
"""

def maximum_fitness(eta):
	if eta < 0.5:
		return 1.0
	else:
		return 1. / (2. * eta)
		
def maximum_abundance(eta):
	if eta < 0.5:
		return 1.0
	else:
		return 1. / float(eta) - 1.
		
fitness_vec = np.vectorize(maximum_fitness)
abundance_vec = np.vectorize(maximum_abundance)

plt.plot(eta_range,fitness_vec(eta_range),lw = 6., label = "Optimal Protocell Composition ($x^*_{FS}$)")
plt.plot(eta_range,abundance_vec(eta_range), lw = 6., label = "Most Abundant Protocell Composition ($\hat{x}_{FS}$)")

plt.xlabel(r"Complementarity Parameter ($\eta)$", fontsize = 24., labelpad = 10.)
plt.ylabel(r"Fraction of Slow Replicators ($x$)", fontsize = 24.)

plt.axvline(x = 0.5, lw = 6., color = "gray", ls = '--', alpha = 0.7)

plt.legend(loc = "lower left", prop = {"size" : 16})

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)



plt.axis([0.0,1.0,0.0,1.04])
plt.tight_layout()


script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)
plt.savefig(protocell_folder + "/Figures/etamodelghostfigure.png")








		
	
	
plt.show()
		