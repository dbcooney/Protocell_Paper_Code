"""
Script used to generate Figure 2.2, illustrating families of densities for various
relative levels of between-protocell competition $\lambda$ for both a case in which
protocell-level competition most favors all-slow compositions (left panel) and a case
in which protocell-level competition most favors an intermediate mix of fast and slow
replicators (right panel).
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


theta = 2.0
s = 0.5

x_step  = 0.0025
x_range = np.arange(x_step, 1.0 +x_step, x_step)


#Formula for density steady states for slow-fast dimorphic protocell model.
	
def steady_state_density(x,lamb,eta,s,theta):
	return (x ** ( (lamb / s) *  (1.0 - eta) - theta - 1.0)) * ((1.0 - x)**(theta - 1.0)) * (np.exp(-(lamb * eta * x)/s))
	
steady_vec = np.vectorize(steady_state_density)

plt.figure(1)

"""
Plotting densities for case in which protocell-level replication rate is maximized by 
all-slow composition.
"""


x_step  = 0.005

eta = .67
density_plot_1 = steady_vec(x_range,4.,eta,s,theta)
density_plot_1 = density_plot_1 / spi.simps(density_plot_1,x_range)
plt.plot(x_range,density_plot_1, lw = 6., color = plt.cm.YlOrRd(0.2), label = r"$\lambda =4 $")

density_plot_2 = steady_vec(x_range,16.,eta,s,theta)
density_plot_2 = density_plot_2/ spi.simps(density_plot_2,x_range)
plt.plot(x_range,density_plot_2, lw = 6., color = plt.cm.YlOrRd(0.4),label = r"$\lambda =16$")

density_plot_3 = steady_vec(x_range,32.,eta,s,theta)
density_plot_3 = density_plot_3 / spi.simps(density_plot_3,x_range)
plt.plot(x_range,density_plot_3, lw = 6., color = plt.cm.YlOrRd(0.6),label = r"$\lambda = 32$")

density_plot_4 = steady_vec(x_range,64.,eta,s,theta)
density_plot_4 = density_plot_4 / spi.simps(density_plot_4,x_range)
plt.plot(x_range,density_plot_4, lw = 6., color = plt.cm.YlOrRd(0.8), label = r"$\lambda =64$")

density_plot_5 = steady_vec(x_range,128.,eta,s,theta)
density_plot_5 = density_plot_5 / spi.simps(density_plot_5,x_range)
plt.plot(x_range,density_plot_5, lw = 6., color = plt.cm.YlOrRd(1.0), label = r"$\lambda =128$")

plt.axis([0.0,1.0,0.0,10.])
plt.legend(loc = "upper right")

plt.xlabel(r"Fraction of Slow Replicators ($x$)", fontsize = 24., labelpad = 20.)
plt.ylabel(r"Probability Density", fontsize = 24.)

plt.axvline(x = 1./(2 * 0.67), lw = 6., ls = '--', color = 'k')
plt.annotate(r"Optimal Composition $x^*_{FS}$",xy = (1./(2 * 0.67) - 0.06,8.),fontsize = 20., rotation = 90.)

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.tight_layout()



script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)
plt.savefig(protocell_folder + "/Figures/etadensitiesghost.png")



plt.figure(2)

"""
Plotting densities for case in which protocell-level replication rate is maximized
by protocell featuring 75 percent slow genes and 25 percent fast genes.
"""


x_step  = 0.0025

eta = .33
density_plot_1 = steady_vec(x_range,2.,eta,s,theta)
density_plot_1 = density_plot_1 / spi.simps(density_plot_1,x_range)
plt.plot(x_range,density_plot_1, lw = 6., color = plt.cm.YlOrRd(0.2), label = r"$\lambda =2$")

density_plot_2 = steady_vec(x_range,4.,eta,s,theta)
density_plot_2 = density_plot_2/ spi.simps(density_plot_2,x_range)
plt.plot(x_range,density_plot_2, lw = 6., color = plt.cm.YlOrRd(0.4), label = r"$\lambda =4$")

density_plot_3 = steady_vec(x_range,8.,eta,s,theta)
density_plot_3 = density_plot_3 / spi.simps(density_plot_3,x_range)
plt.plot(x_range,density_plot_3, lw = 6., color = plt.cm.YlOrRd(0.6), label = r"$\lambda =8$")

density_plot_4 = steady_vec(x_range,16.,eta,s,theta)
density_plot_4 = density_plot_4 / spi.simps(density_plot_4,x_range)
plt.plot(x_range,density_plot_4, lw = 6.,color = plt.cm.YlOrRd(0.8), label = r"$\lambda =16$")

density_plot_5 = steady_vec(x_range,32.,eta,s,theta)
density_plot_5 = density_plot_5 / spi.simps(density_plot_5,x_range)
plt.plot(x_range,density_plot_5, lw = 6., color = plt.cm.YlOrRd(1.), label = r"$\lambda =32$")


plt.axis([0.0,1.0,0.0,10.])
plt.legend(loc = "upper center")

plt.xlabel(r"Fraction of Slow Replicators ($x$)", fontsize = 24., labelpad = 20.)
plt.ylabel(r"Probability Density", fontsize = 24.)


plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.tight_layout()



plt.savefig(protocell_folder + "/Figures/etadensitiesnoghost.png")




plt.show() 