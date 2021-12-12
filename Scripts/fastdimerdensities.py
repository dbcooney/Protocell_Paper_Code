"""
Script used to generate Figure 4.2, illustrating families of steady-state
densities for various relative levels of between-protocell competition $\lambda$ 
for our model of multilevel competition on the fast-dimer edge of the simplex.
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
s = .1

x_step  = 0.0025
x_range = np.arange(x_step, 1.0 +x_step, x_step)



def G_eta(x,eta):
	return 0.5 * (x - 0.5 * eta * (x ** 2.0))
	
def lambda_thresh(eta,theta,s):
	return ( 2.0 * ( s + (1.0/(2.0 +s ))) * theta) / (1.0 - 0.5 * eta)
	
def steady_state_density(x,lamb,eta,s,theta):
	lamb_tilde = (lamb * (2. + s)) / ((s + 1.) ** 2.)
	return (x ** ( lamb_tilde * 0.5 *  (1.0 - 0.5 * eta) - theta - 1.0)) * ((1.0 - x)**(theta - 1.0)) * (np.exp(-(0.25 * lamb_tilde * eta * x)))
	
steady_vec = np.vectorize(steady_state_density)

plt.figure(1)

x_step  = 0.005

"""
Plotting steady state densities for various strengths of between-protocell selection
$\lambda$.
"""

eta = 1.
density_plot_1 = steady_vec(x_range,4.,eta,s,theta)
density_plot_1 = density_plot_1 / spi.simps(density_plot_1,x_range)
plt.plot(x_range,density_plot_1, lw = 6., color = plt.cm.YlOrRd(0.1), label = r"$\lambda =4 $")

density_plot_2 = steady_vec(x_range,8.,eta,s,theta)
density_plot_2 = density_plot_2/ spi.simps(density_plot_2,x_range)
plt.plot(x_range,density_plot_2, lw = 6., color = plt.cm.YlOrRd(0.28),label = r"$\lambda =8$")

density_plot_3 = steady_vec(x_range,16.,eta,s,theta)
density_plot_3 = density_plot_3 / spi.simps(density_plot_3,x_range)
plt.plot(x_range,density_plot_3, lw = 6., color = plt.cm.YlOrRd(0.46),label = r"$\lambda = 16$")

density_plot_4 = steady_vec(x_range,32.,eta,s,theta)
density_plot_4 = density_plot_4 / spi.simps(density_plot_4,x_range)
plt.plot(x_range,density_plot_4, lw = 6., color = plt.cm.YlOrRd(0.64), label = r"$\lambda =32$")

density_plot_5 = steady_vec(x_range,64.,eta,s,theta)
density_plot_5 = density_plot_5 / spi.simps(density_plot_5,x_range)
plt.plot(x_range,density_plot_5, lw = 6., color = plt.cm.YlOrRd(0.82), label = r"$\lambda =64$")

density_plot_6 = steady_vec(x_range,128.,eta,s,theta)
density_plot_6 = density_plot_6 / spi.simps(density_plot_6,x_range)
plt.plot(x_range,density_plot_6, lw = 6., color = plt.cm.YlOrRd(1.0), label = r"$\lambda =128$")

plt.axis([0.0,1.0,0.0,10.])
plt.legend(loc = "upper center")

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.xlabel(r"Fraction of Dimer Replicators $z$", fontsize = 24., labelpad = 10.)
plt.ylabel(r"Probability Density", fontsize = 24.)



plt.tight_layout()

script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)
plt.savefig(protocell_folder + "/Figures/fddensities.png")



plt.show()

