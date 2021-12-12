"""
Script used to generate Figures B.1 and B.2, comparing the analytically calculated 
steady states for the fast-dimer and fast-slow dimorphic models with the states 
achieved after a large number of steps under numerical finite volume simulations of
the same models. 
"""



import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmap
import scipy.integrate as spi
import os


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


"""
This argument toggles between simulating dynamics on the fast-slow or fast-dimer
edges of the simplex.
"""

#edge = "FD"
edge = "FS"


lamb1 = 8.
lamb2 = 16.
lamb3 = 32.
lamb4 =64.
lamb5 = 128.


"""
Introducing protocell-level replication rates for fast-slow and fast-dimer edges.
"""
	
def GFS(x,eta):
	return x - eta * (x ** 2.)
	
def GFD(x,eta):
	return 0.5 * (x - 0.5 * eta * (x ** 2.))
	
def GFS_j(j,N,eta):
	return 0.5 * (1.0 / N) * (2.0 * j + 1.) - (eta / 3.) * (1.0 / (N**2.)) * (3.0 * (j ** 2.) + 3.0 * j + 1.)
	
def GFD_j(j,N,eta):
	return 0.25 * (1. / N) * (2.0 * j + 1.) - (eta / 12.) * (1.0 / (N**2.)) * (3.0 * (j ** 2.) + 3.0 * j + 1.)
	
	
"""
Formulas for analyically calculated steady state densities.
"""	

def fast_dimer_density(x,lamb,eta,s,theta):
	lamb_tilde = (lamb * (2. + s)) / ((s + 1.) ** 2.)
	return (x ** ( lamb_tilde * 0.5 *  (1.0 - 0.5 * eta) - theta - 1.0)) * ((1.0 - x)**(theta - 1.0)) * (np.exp(-(0.25 * lamb_tilde * eta * x)))

def fast_slow_density(x,lamb,eta,s,theta):
	return (x ** ( (lamb / s) *  (1.0 - eta) - theta - 1.0)) * ((1.0 - x)**(theta - 1.0)) * (np.exp(-(lamb * eta * x)/s))
	

def steady_density(x,lamb,eta,s,theta,edge):
	if edge == "FD":
		return fast_dimer_density(x,lamb,eta,s,theta)	
	if edge == "FS":
		return fast_slow_density(x,lamb,eta,s,theta)
	
N = 400
time_step = 0.003
#time_length for FS edge
time_length = 5000
#time_length for FD edge
#time_length = 3000


eta = 0.67
s = 1.
theta = 1.


GFSvec = np.vectorize(GFS)
GFS_jvec = np.vectorize(GFS_j)

GFDvec = np.vectorize(GFD)
GFD_jvec = np.vectorize(GFD_j)

index_holder = np.zeros(N)
for j in range(N):
	index_holder[j] = j
	

"""
Formulating discretized versions of a family of initial densities.
"""

def theta_init(j,N,theta):
	return N ** (1.0 - theta) * (((N - j) ** theta) - ((N - j - 1.0) ** theta) )
	
theta_vec = np.vectorize(theta_init)	
	
f_j = np.ones(N)
f_j = theta_vec(index_holder,N,1.0)



"""
Formulating fluxes between volume boundaries for FS and FD models.
"""

def above_fluxFS(j,N,s):
	return ((j+1.0) / N) * (1.0 - (j+1.0) / N) * s
	
def above_fluxFD(j,N,s):
	return ((j+1.0) / N) * (1.0 - (j+1.0) / N) * (s + 1. / (2. + s))
	
def below_fluxFS(j,N,s):
	return (np.float(j) / N) * (1.0 - (np.float(j)) / N) * s
	
def below_fluxFD(j,N,s):
	return (np.float(j) / N) * (1.0 - (np.float(j)) / N) *  (s + 1. / (2. + s))
	
above_fluxFS_vec = np.vectorize(above_fluxFS)
below_fluxFS_vec = np.vectorize(below_fluxFS)

above_fluxFD_vec = np.vectorize(above_fluxFD)
below_fluxFD_vec = np.vectorize(below_fluxFD)


"""
Characterizing impact of within-protocell dynamics for finite volume discretization.
"""

def within_group(f,N,s,index_holder,edge):
	left_roll = np.roll(f,-1)
	left_roll[-1] = 0.
	right_roll = np.roll(f,1)
	right_roll[0] = 0.
	
	
	if edge == "FS":
		flux_in = above_fluxFS_vec(index_holder,N,s) * left_roll
		flux_out = below_fluxFS_vec(index_holder,N,s) * f
		
	elif edge == "FD":
		flux_in = above_fluxFD_vec(index_holder,N,s) * left_roll
		flux_out = below_fluxFD_vec(index_holder,N,s) * f
		
	return N * (flux_in - flux_out)
	
	
"""
Effect of protocell-level competition on multilevel finite volume dynamics.
"""

def righthand(f,G,N):
	return f * (G - (1.0/N) * np.dot(f,G))
	

peak_holder = [float(np.argmax(f_j))/N]



"""
Setting up finite volume simulations for five different values of relative selection
strength $\lambda$. 
"""
f_j1 = theta_vec(index_holder,N,1.0)
f_j2 = theta_vec(index_holder,N,1.0)
f_j3 = theta_vec(index_holder,N,1.0)
f_j4 = theta_vec(index_holder,N,1.0)
f_j5 = theta_vec(index_holder,N,1.0)



"""
Running the five stes of finite-volume simulations.
"""
for time in range(time_length):

	
	if edge == "FS":
		between_group_effect1 = righthand(f_j1,GFS_jvec(index_holder,N,eta),N)
		between_group_effect2 = righthand(f_j2,GFS_jvec(index_holder,N,eta),N)
		between_group_effect3 = righthand(f_j3,GFS_jvec(index_holder,N,eta),N)
		between_group_effect4 = righthand(f_j4,GFS_jvec(index_holder,N,eta),N)
		between_group_effect5 = righthand(f_j5,GFS_jvec(index_holder,N,eta),N)
	elif edge == "FD":
		between_group_effect1 = righthand(f_j1,GFD_jvec(index_holder,N,eta),N)
		between_group_effect2 = righthand(f_j2,GFD_jvec(index_holder,N,eta),N)
		between_group_effect3 = righthand(f_j3,GFD_jvec(index_holder,N,eta),N)
		between_group_effect4 = righthand(f_j4,GFD_jvec(index_holder,N,eta),N)
		between_group_effect5 = righthand(f_j5,GFD_jvec(index_holder,N,eta),N)
		
	within_group_effect1 = within_group(f_j1,N,s,index_holder,edge)
	within_group_effect2 = within_group(f_j2,N,s,index_holder,edge)
	within_group_effect3 = within_group(f_j3,N,s,index_holder,edge)
	within_group_effect4 = within_group(f_j4,N,s,index_holder,edge)
	within_group_effect5 = within_group(f_j5,N,s,index_holder,edge)
	
	
	righthandside1 = lamb1 * between_group_effect1 + within_group_effect1
	righthandside2 = lamb2 * between_group_effect2 + within_group_effect2
	righthandside3 = lamb3 * between_group_effect3 + within_group_effect3
	righthandside4 = lamb4 * between_group_effect4 + within_group_effect4
	righthandside5 = lamb5 * between_group_effect5 + within_group_effect5
	
	f_j1 = f_j1 + time_step * righthandside1
	f_j2 = f_j2 + time_step * righthandside2
	f_j3 = f_j3 + time_step * righthandside3
	f_j4 = f_j4 + time_step * righthandside4
	f_j5 = f_j5 + time_step * righthandside5
	
	print (1.0 / N) * np.sum(f_j1)
	peak_holder.append(float(np.argmax(f_j1))/N)
	







"""
Plotting the states achieved from the five sets of finite volume simulations after 
the prescribed number of time steps. 
"""

plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j1, color = plt.cm.YlOrRd(0.2), lw = 6., label = r"$\lambda = 8$")
plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j2, color = plt.cm.YlOrRd(0.4), lw = 6., label = r"$\lambda = 16$")
plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j3, color = plt.cm.YlOrRd(0.6), lw = 6., label = r"$\lambda = 32$")
plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j4, color = plt.cm.YlOrRd(0.8), lw = 6., label = r"$\lambda = 64$")
plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j5, color = plt.cm.YlOrRd(1.0), lw = 6., label = r"$\lambda = 128$")


"""
Plotting the analytically calculated steady states for the same values of the 
relative selection strength $\lambda$ as used in the finite volume simulations.
"""

lamb1_steady = steady_density(np.arange(0.5/N,1.0+0.5/N,1.0/N),lamb1,eta,s,1.,edge)
lamb1_norm = lamb1_steady / spi.simps(lamb1_steady,np.arange(0.5/N,1.0+0.5/N,1.0/N))
plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),lamb1_norm, color =plt.cm.YlOrRd(0.2), lw = 8., ls = '--')

lamb2_steady = steady_density(np.arange(0.5/N,1.0+0.5/N,1.0/N),lamb2,eta,s,1.,edge)
lamb2_norm = lamb2_steady / spi.simps(lamb2_steady,np.arange(0.5/N,1.0+0.5/N,1.0/N))
plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),lamb2_norm, color =plt.cm.YlOrRd(0.4), lw = 8., ls = '--')

lamb3_steady = steady_density(np.arange(0.5/N,1.0+0.5/N,1.0/N),lamb3,eta,s,1.,edge)
lamb3_norm = lamb3_steady / spi.simps(lamb3_steady,np.arange(0.5/N,1.0+0.5/N,1.0/N))
plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),lamb3_norm, color =plt.cm.YlOrRd(0.6), lw = 8., ls = '--')

lamb4_steady = steady_density(np.arange(0.5/N,1.0+0.5/N,1.0/N),lamb4,eta,s,1.,edge)
lamb4_norm = lamb4_steady / spi.simps(lamb4_steady,np.arange(0.5/N,1.0+0.5/N,1.0/N))
plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),lamb4_norm, color =plt.cm.YlOrRd(0.8), lw = 8., ls = '--')

lamb5_steady = steady_density(np.arange(0.5/N,1.0+0.5/N,1.0/N),lamb5,eta,s,1.,edge)
lamb5_norm = lamb5_steady / spi.simps(lamb5_steady,np.arange(0.5/N,1.0+0.5/N,1.0/N))
plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),lamb5_norm, color =plt.cm.YlOrRd(1.), lw = 8., ls = '--')


if eta == 0.67:
	plt.axvline(x = 0.75, ls = '--', lw = 6., color = 'k', alpha = 0.8)
	plt.annotate(r"Optimal Composition $x^*_{FS}$", xy = (0.7,8.),rotation = 90., fontsize = 20.)

plt.axis([0.,1.0,0.,10.])

plt.legend(loc = "upper left", prop = {"size":16})


if edge == "FS":	
	plt.xlabel(r"Fraction of Slow Replicators $x$", fontsize = 24., labelpad = 10.)
	plt.ylabel(r"Density of Groups $f(t,x)$", fontsize = 24.)
elif edge == "FD":
	plt.xlabel(r"Fraction of Dimer Replicators $z$", fontsize = 24.,labelpad = 10.)
	plt.ylabel(r"Density of Groups $g(t,z)$", fontsize = 24.)
	
plt.tight_layout()

script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)



if edge == "FS" and eta == 0.67:
	plt.savefig(protocell_folder + "/Figures/fvfsghostdensity.png")
elif edge == "FS" and eta == 0.33:
	plt.savefig(protocell_folder + "/Figures/fvfsnoghostdensity.png")
elif edge == "FD" and eta == 1.:
	plt.savefig(protocell_folder + "/Figures/fvfddensity.png")

plt.show()
