"""
This script is used to generate values that are plotted in the right panel of Figure 5.5.
Here, we calculate numerically the average fraction of slow genes present in the steady 
state densities for the dimorphic fast-dimer competition (with Hölder exponent 
$\theta =1$) for a range of strengths of between-protocell competition $\lambda$.
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi

x_min = 0.
x_max = 1.
x_step = 0.0001
partial_x_range = np.arange(x_min + x_step, x_max + x_step, x_step)
full_x_range = np.arange(x_min, x_max + x_step, x_step)

lamb_min = 0.
lamb_max = 250.
lamb_step = 20.
lamb_range = np.arange(lamb_min,lamb_max + lamb_step, lamb_step)

eta = 1.
s = 1.
theta = 1.

"""
Formulas for the threshold selection strength needed to achieve coexistence, as well as
the formula for steady state densities on the fast-dimer edge and the corresponding
integrand used to compute the mean fraction of slow genes. 
"""

def lamb_tilde(lamb,s):
	b_F = 1 + s
	b_D = (1.0+s) / (2.0 + s)
	lamb_tilde = lamb / (2.0 * (b_F - b_D))
	return lamb_tilde

def steadyFD(lamb, eta, s, z, theta):
	z_exp = 0.5 * lamb_tilde(lamb,s) * (1.0 - 0.5 * eta) - theta - 1.
	if z_exp > -1.:
		return (z ** z_exp) * ((1.0 - z) ** (theta - 1.)) * np.exp(-0.25 *lamb_tilde(lamb,s) * eta * z)
	else:
		return 0.
	
def steadyFD_times_z(lamb, eta, s, z, theta):
	z_exp = 0.5 * lamb_tilde(lamb,s) * (1.0 - 0.5 * eta) - theta
	if z_exp > 0.:
		return (z ** z_exp) * ((1.0 - z) ** (theta - 1.)) * np.exp(-0.25 *lamb_tilde(lamb,s) * eta * z)
	else:
		return 0.
		
steadyFD_vec = np.vectorize(steadyFD)
steadyFD_times_x_vec = np.vectorize(steadyFD_times_z)
	

"""
Calculating the integrals for the mass of groups in the unnormalized steady state density
and the mean fraction of slow replicators for the normalized steady state densities.
"""	
def int_near_zero(lamb,eta,s,z0,theta):
	exp = 0.5 * lamb_tilde(lamb,s) * (1.0 - 0.5 * eta) - theta
	if exp > 0:
		return (1.0 / exp) * (z0 ** (exp))
	else:
		return 0.
	
def mass(lamb,eta,s,z0,theta,x_range):
	near_zero = int_near_zero(lamb,eta,s,x_min,theta)
	rest = spi.simps(steadyFD_vec(lamb,eta,s,partial_x_range,theta),partial_x_range)
	return near_zero + rest
	
def mean(lamb,eta,s,z_0,theta,x_range_1,x_range_2):
	moment_mass = spi.simps(steadyFD_times_x_vec(lamb,eta,s,x_range_2,theta),x_range_2)
	normalizer = mass(lamb,eta,s,x_min,theta,x_range_1)
	if normalizer > 0.:
		return moment_mass / normalizer
	else:
		return 0.
	
print mass(200.,1.,1.,x_min,1.,partial_x_range)
print mean(7200.,1.,1.,x_min,1.,partial_x_range,full_x_range)


"""
Extracting $\lambda$ values used for studying the mean fraction of slow genes under
the trimorphic dynamics.
"""
script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)
file = protocell_folder + "/Simulation_Outputs/meanplottrimorphiceta1.txt"

def process(list):
	list = list.split(',')
	list = [float(a) for a in list]
	return list

f = open(file, 'r+')
lamb_list = f.readline()
#print lamb_list
lamb_list = process(lamb_list)

f.close()


mean_holder = [mean(lamb,eta,s,x_min,theta,partial_x_range,full_x_range) for lamb in lamb_list]
plt.plot(lamb_list,mean_holder, lw = 5., color = 'r', ls = '--')

"""
Saving the numerically computed means of the steady state density for the fast-dimer
competition.
"""
file2 = file_path + "FDmeanploteta1.txt"

f2 = open(file2, 'w+')

f2.write(str(lamb_list)[1:-1])
f2.write('\n')
f2.write(str(mean_holder)[1:-1])
f2.write('\n')
	
f2.close()

plt.show()