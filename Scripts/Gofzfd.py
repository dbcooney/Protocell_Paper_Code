"""
Script used to generate Figure 4.1, illustrating the protocell-level replication
function for dimorphic protocell competition on the fast-dimer edge of the simplex
and various values of the complementarity parameter $\eta$.
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import os

import matplotlib
matplotlib.use('TkAgg')


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

x_step = 0.01
x_range = np.arange(0.,1.+x_step,x_step)


#Defining the protocell-level replication function for fast-dimer competition.

def G(z,eta):
	return 0.5 * (z - 0.5 * eta * (z ** 2.))
	

plt.plot(x_range,G(x_range,0.), lw = 6., label = r"$\eta = 0.0$", color = plt.cm.YlOrRd(0.1))
plt.plot(x_range,G(x_range,0.2), lw = 6., label = r"$\eta = 0.2$", color = plt.cm.YlOrRd(.28))
plt.plot(x_range,G(x_range,0.4), lw = 6., label = r"$\eta = 0.4$", color = plt.cm.YlOrRd(0.46))
plt.plot(x_range,G(x_range,0.6), lw = 6., label = r"$\eta = 0.6$", color = plt.cm.YlOrRd(0.64))
plt.plot(x_range,G(x_range,0.8), lw = 6., label = r"$\eta = 0.8$", color = plt.cm.YlOrRd(0.82))
plt.plot(x_range,G(x_range,1.), lw = 6., label = r"$\eta = 1.0$", color = plt.cm.YlOrRd(1.))

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.xlabel("Fraction of Dimer Replicator $z$", fontsize = 20., labelpad = 10.)
plt.ylabel("Protocell-Level Fitness $G_{FD}(z)$", fontsize = 20.)

plt.legend(loc = "upper left", ncol = 3, prop ={'size':16})

plt.axis([0.0,1.0,0.0,0.625])

plt.tight_layout()


script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)
plt.savefig(protocell_folder + "/Figures/Gofzfdfunction.png")


plt.show()