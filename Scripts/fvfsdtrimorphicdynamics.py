"""
This is our baseline script for running finite volume simulations of the trimorphic
multilevel model for competition between fast, slow, and dimer replicators. 

This script is used to generate Figure 3.2 (illustrating the protocell-level reproduction 
function on the simplex), Figure 5.1 (providing snapshots of the numerical solutions for
trimorphic competition at various points in time), and Figure 5.3 (showing the long-time
states achieved by numerical finite volume simulations for a range of complementarity 
parameters $\eta$ and relative strengths of between-protocell selection $\lambda$.)

This script also generates the data for the average protocell-replication rate as a 
function of time, which is used by the script "avgGintime.py" uses to generate Figure
5.2.
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


s = 1.
bS = 1.0
bF = (1.0 + s)
bD = 1.0 - (1. / (2.0 + s))


"""
Different options for saving images or text data from the simulations: "trajectory"
saves snapshots of the density after a specified time steo, "grouptime" saves a list of 
collective replication rates in the population over the timesteps of the simulation,
and "steady" saves images at the end of long runs (10,000 timesteps) of the simulations.
"""
#quantity = "trajectory"
#quantity = "grouptime"
quantity = "none"
#quantity = "steady"

"""
Setting parameters for multilevel competition model and size of numerical grid.
"""

eta = 0.
N = 100
lamb = 30.


"""
Specifying the time-steppping parameters for the simulations.
"""
timesteps = 10000
time_increments = 0.015



"""
Specifying possible initial distributions.
"""

uniform_dist = np.zeros((N,N))
partial_uniform_dist = np.zeros((N,N))


index_sum = np.zeros((N,N))

for j in range(N):
	for k in range(N):
		index_sum[j,k] = j + k


for j in range(N):
	for k in range(N-j):
		uniform_dist[j,k] = 2.0
		
		if j <= 10 and k <= 10:
			partial_uniform_dist[j,k] = 10.


"""
Discretized protocell-level reproduction function (as derived in Section B.2.1).
"""

def Gjk(j,k,eta,N):

	N = float(N)
	x = float(j) / N
	y = float(k) / N
	
	if j + k < N-1:
		total =  0.5 + 0.5 * x - 0.5 * y - 0.25 * eta - (eta / (24. * (N**2.))) \
		-0.5 * x * eta -0.25 * eta * x * x + 0.5 * y * eta + 0.5 * eta * x * y \
		-0.25 * eta * y * y
		return total
	elif j + k == N-1:
		total = x - eta * x * x + (0.5 / N) - ((7.0 * eta) / (24.0 * N * N)) - ((x * eta) / N )
		return total
	else:
		return None
		
		
Gjk_vec = np.vectorize(Gjk)


def group_birth(j,k,eta,N):
	if j + k < N - 1:
		return Gjk(j,k,eta,N)
	elif j + k == N-1:
		return Gjk(j,k,eta,N)
	else:
		return 0.
		
group_birth_values = np.zeros((N,N))	
	
for j in range(N):
	for k in range(N):
		group_birth_values[j,k] = group_birth(j,k,eta,N)	
		
		
		


spatial_grid = np.zeros((N+1,N+1))
step = 1./float(N)
x_range = np.arange(0.0,1.0 + step, step)
y_range = np.arange(0.0,1.0 + step, step)

def cell_weights(j,k,N):
	if j + k < N - 1:
		return 1.
	elif j + k == N-1:
		return 0.5
	else:
		return 0.
		

		
		
def has_top_right(j,k,N):
	if j+k < N-1:
		return 1.0
	else:
		return 0.0	
		
		
"""
Defining fluxes across the four possible boundaries for a given volume. These fluxes 
correspond to gene-level replication events.
"""		

def left_flux(j,k,N,bS,bF,bD):

	if j + k >= N:
		return 0.

	N = float(N)
	xj = float(j) / N
	yk1 = (k+1.0) / N
	yk = float(k)/ N
	
	coeff = xj * (yk1 - yk)
	inside = (bS - bD) * (1.0 - xj) -  0.5 * (bF - bD) * (yk1 + yk)
	
	
	return coeff * inside
	
	
	
def right_flux(j,k,N,bS,bF,bD):

	if j + k >= N-1:
		return 0.

	N = float(N)
	xj1 = float(j+1.0) / N
	yk1 = (k+1.0) / N
	yk = float(k)/ N
	
	coeff = xj1 * (yk1 - yk)
	inside = (bS - bD) * (1.0 - xj1) -  0.5 * (bF - bD) * (yk1 + yk)
	return -coeff * inside
	
def bottom_flux(j,k,N,bS,bF,bD):

	if j + k >= N:
		return 0.
	
	N = float(N)
	xj = float(j) / N
	xj1 = (j + 1.0) / N
	yk = float(k)/ N
	
	coeff = yk * (xj1 - xj)
	inside = (bF - bD) * (1.0 - yk) - 0.5 * (bS - bD) * (xj1 + xj)
	return coeff * inside
	
	
def top_flux(j,k,N,bS,bF,bD):

	if j + k >= N-1:
		return 0.
	
	N = float(N)
	xj = float(j) / N
	xj1 = (j + 1.0) / N
	yk1 = float(k+1.0)/ N
	
	coeff = yk1 * (xj1 - xj)
	inside = (bF - bD) * (1.0 - yk1) - 0.5 * (bS - bD) * (xj1 + xj)
	return -coeff * inside





cell_weights = np.zeros((N,N))
inv_cell_weights = np.zeros((N,N))

for j in range(N):
	for k in range(N):
		if j + k < N-1:
			cell_weights[j,k] = 1.
			inv_cell_weights[j,k] = 1.
		elif j + k == N-1:
			cell_weights[j,k] = 0.5
			inv_cell_weights[j,k] = 2.
		


def group_righthand(state,group_birth_values,eta,N,cell_weights):
	weighted = np.multiply(state,group_birth_values)
	increase = weighted
	decrease = (1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights)) * state
	return increase-decrease


"""
Calculating fluxes for the boundaries of each volume in our grid.
"""	
left_flux_array = np.zeros((N,N))
right_flux_array = np.zeros((N,N))
bottom_flux_array = np.zeros((N,N))
top_flux_array = np.zeros((N,N))

for j in range(N):
	for k in range(N):
		left_flux_array[j,k] = left_flux(j,k,N,bS,bF,bD)
		right_flux_array[j,k] = right_flux(j,k,N,bS,bF,bD)
		bottom_flux_array[j,k] = bottom_flux(j,k,N,bS,bF,bD)
		top_flux_array[j,k] = top_flux(j,k,N,bS,bF,bD)
	
		
"""
Choosing initial condition.
"""		
	
state = uniform_dist
#state = partial_uniform_dist


cell_indices = np.zeros((N,N,2))
for j in range(N):
	for k in range(N):
		cell_indices[j,k,0] = j
		cell_indices[j,k,1] = k



"""
Calculating impact of gene-level competition dynamics on finite volume representation
of the densities. 
"""	

def within_group(state,left_flux_array,right_flux_array,bottom_flux_array,top_flux_array,cell_weights):
	
	up_roll = np.roll(state,1,axis = 1)
	up_roll[:,0] = 0.
	
	down_roll = np.roll(state,-1,axis = 1)
	down_roll[:,N-1] = 0.
	
	left_roll = np.roll(state,-1,axis = 0)
	left_roll[N-1,0] = 0.
	
	right_roll = np.roll(state,1,axis = 0)
	right_roll[0,:] = 0.
	
	top_flux_up = np.ones(N)
	bottom_flux_up = np.ones(N)

	
	right_flux_from_right = np.where(right_flux_array > 0.,1.0,0.0)
	right_flux_from_left = np.where(right_flux_array < 0.,1.0,0.0)
	
	left_flux_from_right = np.where(left_flux_array < 0.,1.0,0.0)
	left_flux_from_left = np.where(left_flux_array > 0.,1.0,0.0)
	
	top_contribution = top_flux_array * state
	bottom_contribution = bottom_flux_array * up_roll
	
	right_contribution = right_flux_array * right_flux_from_right * left_roll + right_flux_array * right_flux_from_left * state
	left_contribution = left_flux_array * left_flux_from_right * state + left_flux_array * left_flux_from_left * right_roll
	
	sum = bottom_contribution + top_contribution + left_contribution + right_contribution
	
	
	
	weighted_sum = np.multiply(sum,inv_cell_weights)
	
	return (N**2.) * weighted_sum 

		


"""
Heatmap plotting protocell-level reproduction function $G(x,y,z)$ on the simplex.
"""

	
G_values = np.zeros((N+1,N+1))

for j in range(N+1):
	for k in range(N+1):
		G_values[j,k] = Gjk(j,k,eta,N)
		


cmap = plt.get_cmap('jet')
#cmap = plt.jet()
cmap.set_bad('w',None)

G_values = np.fliplr(G_values)
G_values = np.transpose(G_values)

plt.imshow(G_values)
plt.colorbar(pad = 0.02)

x_ticks = [0.0,0.2,0.4,0.6,0.8,1.0]
y_ticks = [1.0,0.8,0.6,0.4,0.2,0.0]
plt.xticks(range(0,N+1,20),x_ticks)
plt.yticks(range(0,N+1,20),y_ticks)

plt.xticks(fontsize = 14., rotation = 0)
plt.yticks(fontsize = 14.)





plt.tick_params(top = False, right = False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)


script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)


if eta == 0.7:
	plt.savefig(protocell_folder + "/Figures/Gofxyeta0p7first.png",transparent = True)
	plt.ylabel(r"Fraction of Fast Replicators $y$", fontsize = 20.)
	plt.xlabel(r"Fraction of Slow Replicators $x$", fontsize = 20.,labelpad = 10.)
	plt.tight_layout()
	plt.savefig(protocell_folder + "/Figures/Gofxyeta0p7second.png",transparent = True)
elif eta == 1.0:
	plt.ylabel(r"Fraction of Fast Replicators $y$", fontsize = 20.)
	plt.savefig(protocell_folder + "/Figures/Gofxyeta1first.png",transparent = True)
	plt.xlabel(r"Fraction of Slow Replicators $x$", fontsize = 20.,labelpad = 10.)
	plt.tight_layout()
	plt.savefig(protocell_folder + "/Figures/Gofxyeta1second.png",transparent = True)
	
elif eta == 0.9:
	plt.savefig(protocell_folder + "/Figures/Gofxyeta0p9first.png",transparent = True)
	plt.xlabel(r"Fraction of Slow Replicators $x$", fontsize = 20.,labelpad = 10.)
	plt.ylabel(r"Fraction of Fast Replicators $y$", fontsize = 20.)
	plt.tight_layout()
	plt.savefig(protocell_folder + "/Figures/Gofxyeta0p9second.png",transparent = True)
	
	
elif eta == 0.0:
	plt.savefig(protocell_folder + "/Figures/Gofxyeta0first.png",transparent = True)
	plt.xlabel(r"Fraction of Slow Replicators $x$", fontsize = 20.,labelpad = 10.)
	plt.ylabel(r"Fraction of Fast Replicators $y$", fontsize = 20.)
	plt.tight_layout()
	plt.savefig(protocell_folder + "/Figures/Gofxyeta0second.png",transparent = True)

plt.figure(2)


avgGholder = []

"""
Running the finite volume simulation for the trimorphic multilevel competition.
"""

for time in range(timesteps):
	old_state = state 
	 
	within = within_group(state,left_flux_array,right_flux_array,bottom_flux_array,top_flux_array,inv_cell_weights)
	between = group_righthand(state,group_birth_values,eta,N,cell_weights)
	righthand =  within + lamb * between 
	state = state + time_increments * righthand
	
	if time % 40 == 0:
		print np.sum(np.multiply(state,cell_weights))

	
	print np.sum(np.multiply(state,cell_weights))
	print np.sum(np.multiply(state - old_state,cell_weights))
	print np.multiply(state - old_state,cell_weights)
	weighted = np.multiply(state,group_birth_values)
	avgGholder.append((1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights)))
	
	print (1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights)) 

weighted = np.multiply(state,group_birth_values)
print  (1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights))
print 0.5 * (1.0 - 0.5 * eta)
print 1 - eta

"""
Making heat map of the final state achieved as a numerical solution for the multilevel
competition.
"""

state = np.where(index_sum >= N,np.nan,state)
state = np.fliplr(state)
state = np.transpose(state)




cmap = plt.get_cmap('jet')
#cmap = plt.jet()
cmap.set_bad('w',None)



plt.imshow(state,cmap =cmap)
plt.colorbar()

x_ticks = [0.0,0.2,0.4,0.6,0.8,1.0]
y_ticks = [1.0,0.8,0.6,0.4,0.2,0.0]
plt.xticks(range(0,N+1,20),x_ticks)
plt.yticks(range(0,N+1,20),y_ticks)
print np.arange(0.,1.0 + .1,.1)

plt.xticks(fontsize = 14., rotation = 0)
plt.yticks(fontsize = 14.)

#plt.tight_layout()

plt.tick_params(top = False, right = False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)



"""
Plotting the long-time steady state and saving the resulting image.
"""

if quantity == "steady":
	if lamb == 10. and eta == 1.0:
		plt.ylabel(r"Fraction of Fast Replicators $y$", fontsize = 24.)
		plt.title(r"$\eta = 1.0$, $\lambda = 10$", fontsize = 24.)
		#plt.tight_layout()
		plt.savefig(protocell_folder + "/Figures/fvsteadylamb10eta1.png",transparent = True,bbox_inches='tight',pad = 0)
	elif lamb == 30. and eta == 1.0:
		plt.title(r"$\eta = 1.0$, $\lambda = 30$", fontsize = 24.)
		plt.ylabel(r"Fraction of Fast Replicators $y$", fontsize = 24.)
		plt.savefig(protocell_folder + "/Figures/fvsteadylamb30eta1.png",transparent = True,bbox_inches='tight',pad = 0)
		#plt.tight_layout()
	elif lamb == 50. and eta == 1.0:
		plt.title(r"$\eta = 1.0$, $\lambda = 50$", fontsize = 24.)
		plt.ylabel(r"Fraction of Fast Replicators $y$", fontsize = 24.)
		#plt.tight_layout()
		plt.savefig(protocell_folder + "/Figures/fvsteadylamb50eta1.png",transparent = True,bbox_inches='tight',pad = 0)
	elif lamb == 100. and eta == 1.0:
		plt.title(r"$\eta = 1.0$, $\lambda = 100$", fontsize = 24.)
		plt.xlabel(r"Fraction of Slow Replicators $x$", fontsize = 24.,labelpad = 10.)
		#plt.tight_layout()
		plt.ylabel(r"Fraction of Fast Replicators $y$", fontsize = 24.)
		plt.savefig(protocell_folder + "/Figures/fvsteadylamb100eta1.png",transparent = True,bbox_inches='tight',pad = 0)
	elif lamb == 10. and eta == 0.7:
		plt.title(r"$\eta = 0.7$, $\lambda = 10$", fontsize = 24.)
		plt.savefig(protocell_folder + "/Figures/fvsteadylamb10eta0p7.png",transparent = True,bbox_inches='tight',pad = 0)	
	elif lamb == 30. and eta == 0.7:
		plt.title(r"$\eta = 0.7$, $\lambda = 30$", fontsize = 24.)
		plt.savefig(protocell_folder + "/Figures/fvsteadylamb30eta0p7.png",transparent = True,bbox_inches='tight',pad = 0)
	elif lamb == 50. and eta == 0.7:
		plt.title(r"$\eta = 0.7$, $\lambda = 50$", fontsize = 24.)
		plt.savefig(protocell_folder + "/Figures/fvsteadylamb50eta0p7.png",transparent = True,bbox_inches='tight',pad = 0)
	elif lamb == 100. and eta == 0.7:
		plt.title(r"$\eta = 0.7$, $\lambda = 100$", fontsize = 24.)
		plt.xlabel(r"Fraction of Slow Replicators $x$", fontsize = 24.,labelpad = 10.)
		plt.savefig(protocell_folder + "/Figures/fvsteadylamb100eta0p7.png",transparent = True,bbox_inches='tight',pad = 0)
	elif lamb == 10. and eta == 0.9:
		plt.title(r"$\eta = 0.9$, $\lambda = 10$", fontsize = 24.)
		plt.savefig(protocell_folder + "/Figures/fvsteadylamb10eta0p9.png",transparent = True,bbox_inches='tight',pad = 0)
	elif lamb == 50. and eta == 0.9:
		plt.title(r"$\eta = 0.9$, $\lambda = 50$", fontsize = 24.)
		plt.savefig(protocell_folder + "/Figures/fvsteadylamb50eta0p9.png",transparent = True,bbox_inches='tight',pad = 0)
	elif lamb == 100. and eta == 0.9:
		plt.title(r"$\eta = 0.9$, $\lambda = 100$", fontsize = 24.)
		plt.xlabel(r"Fraction of Slow Replicators $x$", fontsize = 24.,labelpad = 10.)
		plt.savefig(protocell_folder + "/Figuresfvsteadylamb100eta0p9.png",transparent = True,bbox_inches='tight',pad = 0)	
	elif lamb == 30. and eta == 0.9:
		plt.title(r"$\eta = 0.9$, $\lambda = 30$", fontsize = 24.)
		plt.savefig(protocell_folder + "/Figures/fvsteadylamb30eta0p9.png",transparent = True,bbox_inches='tight',pad = 0)


"""
Plotting numerically computed solutions after short number of time steps, and saving the
images for Figure 5.1.
"""

elif quantity == "trajectory":

	plt.xlabel(r"Fraction of Slow Replicators $x$", fontsize = 24.,labelpad = 10.)
	plt.ylabel(r"Fraction of Fast Replicators $y$", fontsize = 24.)
	
	if timesteps == 10:
		plt.title(r"$10$ timesteps", fontsize = 24.)
		plt.savefig(protocell_folder + "/Figures/fvtrajectorytime10.png",transparent = True,bbox_inches='tight',pad = 0)
	if timesteps == 20:
		plt.title(r"$20$ timesteps", fontsize = 24.)
		plt.savefig(protocell_folder + "/Figures/fvtrajectorytime20.png",transparent = True,bbox_inches='tight',pad = 0)
	if timesteps == 50:
		plt.title(r"$50$ timesteps", fontsize = 24.)
		plt.savefig(protocell_folder + "/Figures/fvtrajectorytime50.png",transparent = True,bbox_inches='tight',pad = 0)
	if timesteps == 100:
		plt.title(r"$100$ timesteps", fontsize = 24.)
		plt.savefig(protocell_folder + "/Figures/fvtrajectorytime100.png",transparent = True,bbox_inches='tight',pad = 0)
	if timesteps == 150:
		plt.title(r"$150$ timesteps", fontsize = 24.)
		plt.savefig(protocell_folder + "/Figures/fvtrajectorytime150.png",transparent = True,bbox_inches='tight',pad = 0)
	if timesteps == 200:
		plt.title(r"$200$ timesteps", fontsize = 24.)
		plt.savefig(protocell_folder + "/Figures/fvtrajectorytime200.png",transparent = True,bbox_inches='tight',pad = 0)
	
	


"""
Saving list of average protocell-level replication rates over the numerical solution for 
the density at each step in time. This data is then used to generate Figure 5.2.
"""


if quantity == "grouptime":
	file1 = protocell_folder + "/Simulation_Outputs/group_avg_time.txt"
	f1 = open(file1,'a+')
	if lamb == 100 and eta == 1.:
		f1.write("timesteps")
		f1.write('\n')
		f1.write(str(timesteps))
		f1.write('\n')
	f1.write("lambda")
	f1.write('\n')
	f1.write(str(lamb))
	f1.write('\n')
	f1.write("eta")
	f1.write('\n')
	f1.write(str(eta))
	f1.write('\n')
	f1.write("average of G")
	f1.write('\n')
	f1.write(str(avgGholder)[1:-1])
	f1.write('\n')
	f1.close()

	

plt.show()
		
		