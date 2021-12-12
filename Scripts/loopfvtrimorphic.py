"""
This is a script for running finite volume simulations of the trimorphic
multilevel model (for competition between fast, slow, and dimer replicators) over a range
of values for the relative selection strength $\lambda$. 


This script generates the data used to produce Figures 5.4 and 5.5, which respectively
illustrate the collective outcome (in terms of average protocell-level replication rate)
and the mean and modal fraction of slow genes after a large number of time steps for
our finite volume simulation, plotted as a function of the relative strength of
protocell-level replication $\lambda$.
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
#bD = 1.


"""
Specifying options of which data we would like to save during a given simulation run. The
options are to save the protocell-level fitness, the fraction of slow replicators in the
modal composition, and the average fraction of slow genes for the long-run state achieved
by the numerical simulation.
"""
quantity = "fitness"
#quantity = "peak"
#quantity = "mean"


eta = 1.0


N = 100

lamb_min = 0.
lamb_max = 250.
lamb_step = 0.5
lamb_range = np.arange(lamb_min,lamb_max + lamb_step, lamb_step)

lamb = 5.

timesteps = 5000
time_increments = 0.015

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

spatial_grid = np.zeros((N+1,N+1))
step = 1./float(N)
x_range = np.arange(0.0,1.0 + step, step)
y_range = np.arange(0.0,1.0 + step, step)


"""
Defining protocell-level replication rates (as derived in Section B.2.1.)
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
		return 0
		
		
Gjk_vec = np.vectorize(Gjk)


def cell_weights(j,k,N):
	if j + k < N - 1:
		return 1.
	elif j + k == N-1:
		return 0.5
	else:
		return 0.
		

		


def group_birth(j,k,eta,N):
	if j + k < N - 1:
		return Gjk(j,k,eta,N)
	elif j + k == N-1:
		return Gjk(j,k,eta,N)
	else:
		return 0.
		
group_birth_values = np.zeros((N,N))
slow_values = np.zeros((N,N))
fast_values = np.zeros((N,N))	
	
for j in range(N):
	for k in range(N):
		group_birth_values[j,k] = group_birth(j,k,eta,N)
		jx = np.float(j)
		ky = np.float(k)
		if j + k < N - 1:
			slow_values[j,k] = jx / N + 1.0 / (2.0 * N)
			fast_values[j,k] = ky / N + 1.0 / (2.0 * N)
		elif j+k == N - 1:
			slow_values[j,k] = jx / N + 1.0 / (3.0 * N)
			fast_values[j,k] = (N - jx) / N - 1.0 / (3.0 * N)
			
			

"""
Defining effect of between-protocell competition for finite volume simulation.
"""

def group_righthand(state,group_birth_values,eta,N,cell_weights):
	weighted = np.multiply(state,group_birth_values)
	increase = weighted
	
	decrease = (1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights)) * state
	return increase-decrease
		
		
def has_top_right(j,k,N):
	if j+k < N-1:
		return 1.0
	else:
		return 0.0	
		

"""
Specifying fluxes across volume boundaries (corresponding to gene-level dynamics).
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




print group_birth_values

print np.size(uniform_dist)
print np.size(group_birth_values)

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
		


		
		
"""
Calculating the fraction of slow and fast replicators across population of protocells.
"""		
def slow_replicator(state,N,cell_weights):
	weighted = np.multiply(state,slow_values)
	return (1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights))
	
		
def fast_replicator(state,N,cell_weights):
	weighted = np.multiply(state,fast_values)
	return (1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights))

"""
Setting initial distribution.
"""	
state = uniform_dist


cell_indices = np.zeros((N,N,2))
for j in range(N):
	for k in range(N):
		cell_indices[j,k,0] = j
		cell_indices[j,k,1] = k


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
Defining gene-level dynamics for finite volume simulation.
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

	
G_values = np.zeros((N+1,N+1))

for j in range(N+1):
	for k in range(N+1):
		G_values[j,k] = Gjk(j,k,1.0,N)
		
print G_values
G_values = np.flipud(G_values)
print G_values

plt.imshow(G_values)
plt.colorbar()

plt.figure(2)


"""
Calculating the comparable steady state collective outcomes for dimorphic competition
on the fast-dimer edge of the simplex (as derived in Section 4.1).
"""

def Group_FD(lamb,theta,eta,s):
	G1 = 0.5 * (1.0 - 0.5 * eta)
	indv = s + ( 1. / (2.0 + s))
	
	if lamb * G1 <= indv * theta:
		return 0
	
	else:
		return G1 - ((indv * theta) / lamb)

		
GroupFDvec = np.vectorize(Group_FD)

"""
Calculating fraction of slow genes on a given volume.
"""
def slow_material(x,y):
	return x + 0.5 * (1. - x - y)


lamb_holder_list = []
Group_holder = []
Group_FD_holder = []
peak_holder = []
average_comp_holder = []

x_mode = []
y_mode = []


"""
Running finite volume simulations for a range of levels of between-protocell competition.
"""

for lamb_holder in lamb_range:
	lamb = lamb_holder
	state = uniform_dist
	for time in range(timesteps):
	 
		within = within_group(state,left_flux_array,right_flux_array,bottom_flux_array,top_flux_array,inv_cell_weights)
		between = group_righthand(state,group_birth_values,eta,N,cell_weights)
		righthand =  within + lamb * between 
		state = state + time_increments * righthand
	
		weighted = np.multiply(state,group_birth_values)
		
		
	lamb_holder_list.append(lamb)
	Group_holder.append((1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights)))
	Group_FD_holder.append(Group_FD(lamb,1.0,eta,s))
	

	modal = np.unravel_index(np.argmax(state),np.shape(state))
	print modal
	modal = list(modal)
	peak_holder.append(slow_material(np.float(modal[0]) / N, np.float(modal[1])/N))
	
	x_mode.append(np.float(modal[0]) / N)
	y_mode.append(np.float(modal[1]) / N)
	
	average_comp = slow_material(slow_replicator(state,N,cell_weights),fast_replicator(state,N,cell_weights))
	average_comp_holder.append(average_comp)
	
	print lamb



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

plt.xlabel(r"Fraction of Slow Replicators $x$", fontsize = 20.,labelpad = 10.)
plt.ylabel(r"Fraction of Fast Replicators $y$", fontsize = 20.)

plt.tight_layout()




plt.figure(3)

print Group_holder
plt.plot(lamb_range,Group_holder, lw = 5.)
print Group_FD(200.,1.,eta,s)
print GroupFDvec(lamb_range,1.0,eta,s)
plt.plot(lamb_range,Group_FD_holder, color = 'g', lw = 5., ls = "--")

plt.axhline(y = 0.5 * (1.0 - 0.5 * eta), lw = 5.,ls = '--', color = 'k', alpha = 0.9)

plt.axis([0,lamb_max + 1.,0.,0.4 ])

plt.xlabel(r"Cell-Level Selection Strength ($\lambda$)", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Average Cell-Level Reproduction Rate ($\langle G \rangle_f$)", fontsize = 20.)

plt.figure(4)

print peak_holder

plt.plot(lamb_range,peak_holder, lw = 5.)
plt.plot(lamb_range,average_comp_holder, lw = 5.)

plt.xlabel(r"Cell-Level Selection Strength ($\lambda$)", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Fraction of Slow Replicators ($x + \frac{z}{2}$)", fontsize = 20.)


"""
Saving data for producing Figures 5.4 and 5.5.
"""

script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)
file_path = protocell_folder




if quantity == "fitness":
	if eta == 1.0:
		file = file_path + "/Simulation_Outputs/Gplottrimorphiceta1.txt"
	elif eta == 0.9:
		file = file_path + "/Simulation_Outputs/Gplottrimorphiceta0p9.txt"
	elif eta == 0.8:
		file = file_path + "/Simulation_Outputs/Gplottrimorphiceta0p8.txt"
	elif eta == 0.7:
		file = file_path + "/Simulation_Outputs/Gplottrimorphiceta0p7.txt"
	f = open(file,'w+')
	f.write(str(lamb_holder_list)[1:-1])
	f.write('\n')
	f.write(str(Group_holder)[1:-1])
	f.write('\n')
	f.write(str(Group_FD_holder)[1:-1])
	f.write('\n')
	f.close()
	
elif quantity == "peak":
	if eta == 1.0:
		file = file_path + "/Simulation_Outputs/peakplottrimorphiceta1.txt"
	
	f = open(file, 'w+')
	f.write(str(lamb_holder_list)[1:-1])
	f.write('\n')
	f.write(str(peak_holder)[1:-1])
	f.write('\n')
	f.write(str(x_mode)[1:-1])
	f.write('\n')
	f.write(str(y_mode)[1:-1])
	f.write('\n')
	
	f.close()
	
elif quantity == "mean":
	if eta == 1.0:
		file = file_path + "/Simulation_Outputs/
		meanplottrimorphiceta1.txt"
	
	f = open(file, 'w+')
	f.write(str(lamb_holder_list)[1:-1])
	f.write('\n')
	f.write(str(peak_holder)[1:-1])
	f.write('\n')
	
	f.close()

plt.show()
		
		