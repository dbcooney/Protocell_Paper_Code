"""
Script for finite volume numerical simulations for dimorphic multilevel competition. 
Used to generate Figure 4.3, which illustrates time-dependent numerical solutions to
the dynamics on the fast-slow edge of the simplex (left panel) and on the fast-dimer
edge of the simplex (right panel), starting from a uniform initial density.
"""




import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmap
import os


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)



"""
Options to run dynamics for fast-slow or fast-dimer edges of the simplex.
"""
#edge = "FS"
edge = "FD"



"""
Formulating protocell-level replication rates for fast-slow and fast-dimer edges of the
simplex.
"""	
def GFS(x,eta):
	return x - eta * (x ** 2.)
	
def GFD(x,eta):
	return 0.5 * (x - 0.5 * eta * (x ** 2.))
	
def GFS_j(j,N,eta):
	return 0.5 * (1.0 / N) * (2.0 * j + 1.) - (eta / 3.) * (1.0 / (N**2.)) * (3.0 * (j ** 2.) + 3.0 * j + 1.)
	
def GFD_j(j,N,eta):
	return 0.25 * (1. / N) * (2.0 * j + 1.) - (eta / 12.) * (1.0 / (N**2.)) * (3.0 * (j ** 2.) + 3.0 * j + 1.)


GFSvec = np.vectorize(GFS)
GFS_jvec = np.vectorize(GFS_j)

GFDvec = np.vectorize(GFD)
GFD_jvec = np.vectorize(GFD_j)

	
N = 800
time_step = 0.003
#time_length for FS edge
#time_length = 1800
#time_length for FD edge
time_length = 3000


eta = 1.
s = 1.
lamb = 20.



"""
Different possible families of initial densities.
"""

def theta_init(j,N,theta):
	return N ** (1.0 - theta) * (((N - j) ** theta) - ((N - j - 1.0) ** theta) )
	
theta_vec = np.vectorize(theta_init)


index_holder = np.zeros(N)
for j in range(N):
	index_holder[j] = j
	
f_j = np.ones(N)
f_j = theta_vec(index_holder,N,1.0)
#f_j = np.zeros(N)


"""
Formulating within-protocell dynamics as fluxes across volume boundaries.
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
Computing effects of within-protocell competition.
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
Computing effects of between-protocell competition.
"""

def righthand(f,G,N):
	return f * (G - (1.0/N) * np.dot(f,G))

peak_holder = [float(np.argmax(f_j))/N]



"""
Running finite volume simulations for multilevel dynamics.
"""

for time in range(time_length):

	if time % 100 == 0:
		plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j, color = cmap.jet((np.float(time) / time_length)**.25), lw = 3.)
	elif time in [20,40,60,80,100]:
		plt.plot(np.arange(0.5/N,1.0+0.5/N,1.0/N),f_j, color = cmap.jet((np.float(time) / time_length)**.25), lw = 3.)

	if edge == "FS":
		between_group_effect = righthand(f_j,GFS_jvec(index_holder,N,eta),N)
	elif edge == "FD":
		between_group_effect = righthand(f_j,GFD_jvec(index_holder,N,eta),N)
		
	within_group_effect = within_group(f_j,N,s,index_holder,edge)
	righthandside = lamb * between_group_effect + within_group_effect
	f_j = f_j + time_step * righthandside
	
	print (1.0 / N) * np.sum(f_j)
	peak_holder.append(float(np.argmax(f_j))/N)
	



if edge == "FS":	
	plt.xlabel(r"Fraction of Slow Replicators ($x$)", fontsize = 20.)
	plt.ylabel(r"Density of Groups ($f(t,x)$)", fontsize = 20.)
elif edge == "FD":
	plt.xlabel(r"Fraction of Dimer Replicators ($z$)", fontsize = 20.)
	plt.ylabel(r"Density of Groups ($g(t,z)$)", fontsize = 20.)
	

plt.tight_layout()



script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)

if edge == "FS":
	plt.savefig(protocell_folder + "/Figures/FSsampletrajectory.png")
elif edge == 'FD':
	plt.savefig(protocell_folder + "/Figures/FDsampletrajectory.png")





plt.show()
