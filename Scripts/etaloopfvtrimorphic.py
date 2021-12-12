"""
Script used to generate Figure 5.6, illustrating the average protocell-level fitness
in the states achieved after a long-time run of the finite volume simuation, plotted as
a function of the complementarity parameter $\eta$. The figure also provides a
comparison with the analytically calculated collective outcomes achieved on the 
fast-slow and fast-dimer edges of the simplex.
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

eta = 0.8
N = 100

"""
Range of complementarity parameters used for simulations.
"""
eta_min = 0.0
eta_max = 1.0
eta_step = 0.05
eta_range = np.arange(eta_min,eta_max + eta_step, eta_step)


lamb = 75.


timesteps = 5000
time_increments = 0.015

uniform_dist = np.zeros((N,N))
partial_uniform_dist = np.zeros((N,N))


index_sum = np.zeros((N,N))


for j in range(N):
	for k in range(N):
		index_sum[j,k] = j + k


"""
Listing possible initial conditions.
"""

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
Defining protocell-level replication rates.
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
		#return 0.5* Gjk(j,k,eta,N)
		return Gjk(j,k,eta,N)
	else:
		return 0.
		
group_birth_values = np.zeros((N,N))
slow_values = np.zeros((N,N))
fast_values = np.zeros((N,N))	
	
for j in range(N):
	for k in range(N):
		group_birth_values[j,k] = group_birth(j,k,eta,N)
		slow_values[j,k] = np.float(j) / N	
		fast_values[j,k] = np.float(k) / N
	
	
		
		
def has_top_right(j,k,N):
	if j+k < N-1:
		return 1.0
	else:
		return 0.0	
	

"""
Defining fluxes across volume boundaries (corresponding to gene-level dynamics).
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
			

def group_righthand(state,group_birth_values,eta,N,cell_weights):
	weighted = np.multiply(state,group_birth_values)

	increase = weighted
	decrease = (1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights)) * state

	return increase-decrease
		
def slow_replicator(state,N,cell_weights):
	weighted = np.multiply(state,slow_values)
	return (1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights))
	
		
def fast_replicator(state,N,cell_weights):
	weighted = np.multiply(state,fast_values)
	return (1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights))
	
state = uniform_dist

print np.sum(state)
print np.sum(np.multiply(state,cell_weights))

cell_indices = np.zeros((N,N,2))
for j in range(N):
	for k in range(N):
		cell_indices[j,k,0] = j
		cell_indices[j,k,1] = k



"""
Describing effects of gene-level dynamics.
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
Calculating collective outcomes for steady states from dimorphic multilevel competition
on fast-slow and fast-dimer edges of the simplex.
"""

def Group_FD(lamb,theta,eta,s):
	G1 = 0.5 * (1.0 - 0.5 * eta)
	indv = s + ( 1. / (2.0 + s))
	
	if lamb * G1 <= indv * theta:
		return 0
	else:
		return G1 - ((indv * theta) / lamb)
	
		
def Group_FS(lamb,theta,eta,s):
	G1 = 1.0 - eta
	indv = s
	if lamb * G1 <= indv * theta:
		return 0.
	else:
		return G1 - ((indv * theta) / lamb)
		
def Group_max(eta):
	if eta < 0.5:
		return 1.0 - eta
	else:
		return 1. / (4.0 * eta)
		
GroupFDvec = np.vectorize(Group_FD)

def slow_material(x,y):
	return x + 0.5 * (1. - x - y)


eta_holder_list = []
Group_holder = []
Group_FD_holder = []
Group_FS_holder = []
peak_holder = []
average_comp_holder = []


"""
Running finite volume simulations for trimorphic dynamics across the range of possible
complementarity parameters $\eta$.
"""

for eta_holder in eta_range:
	eta = eta_holder
	state = uniform_dist
	
	group_birth_values = np.zeros((N,N))
	slow_values = np.zeros((N,N))
	fast_values = np.zeros((N,N))	
	
	for j in range(N):
		for k in range(N):
			group_birth_values[j,k] = group_birth(j,k,eta,N)
			slow_values[j,k] = np.float(j) / N	
			fast_values[j,k] = np.float(k) / N
	
	
	for time in range(timesteps):
	 
		within = within_group(state,left_flux_array,right_flux_array,bottom_flux_array,top_flux_array,inv_cell_weights)
		between = group_righthand(state,group_birth_values,eta,N,cell_weights)
		righthand =  within + lamb * between 
		state = state + time_increments * righthand
	
	
		weighted = np.multiply(state,group_birth_values)
		
		
	
	weighted = np.multiply(state,group_birth_values)
		
	eta_holder_list.append(eta)
	Group_holder.append((1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights)))
	Group_FD_holder.append(Group_FD(lamb,1.0,eta,s))
	Group_FS_holder.append(Group_FS(lamb,1.0,eta,s))
	
	
	modal = np.unravel_index(np.argmax(state),np.shape(state))
	print modal
	modal = list(modal)
	peak_holder.append(slow_material(np.float(modal[0]) / N, np.float(modal[1])/N))
	
	average_comp = slow_material(slow_replicator(state,N,cell_weights),fast_replicator(state,N,cell_weights))
	average_comp_holder.append(average_comp)
	
	print lamb


print "Cheesecake"
weighted = np.multiply(state,group_birth_values)
print  (1.0 / (float(N) ** 2.)) * np.sum(np.multiply(weighted,cell_weights))
print 0.5 * (1.0 - 0.5 * eta)
print 1 - eta

state = np.where(index_sum >= N,np.nan,state)
state = np.fliplr(state)
state = np.transpose(state)



print index_sum
print state[N-2,N-2]
print "test"


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

"""
Plotting comparison between collective outcome at steady state for trimorphic
numerical simulations and analytically calculated collective fitness for steady states
on the fast-slow and fast-dimer edges of the simplex.
"""

plt.plot(eta_range,Group_holder, lw = 7., label = r"$\langle G(x,y) \rangle_{u(x,y)}$ (Trimorphic)")


Group_max_holder = [Group_max(eta) for eta in eta_range]

if lamb == 10.:
	plt.plot(eta_range,Group_FD_holder, color = 'r', lw = 7., ls = "-.", label = r"$\langle G_{FD} \rangle_{g_{FD}(z)}$, (dimorphic, FD edge)")
	

G_alldimer = [0.5 * (1.0 - 0.5 * eta) for eta in eta_range]
G_allslow = [1.0 - eta for eta in eta_range]
plt.plot(eta_range,G_alldimer, color = 'r', lw = 7., ls = '--', alpha = 1., label = "$G(0,0)$: All-Dimer Fitness")

if lamb == 10.:
	plt.plot(eta_range,Group_FS_holder, color = 'k', lw = 7., ls = "-.", label = r"$\langle G_{FD} \rangle_{g_{FS}(x)}$, (dimorphic, FS edge)")

plt.plot(eta_range,G_allslow, color = 'k', lw = 7., ls = '--', alpha = 1., label = "$G(1,0)$ : All-Slow Fitness")


plt.xlabel(r"Complementarity Parameter ($\eta$)", fontsize = 20., labelpad = 10.)
plt.ylabel(r"Average Protocell-Level Fitness ($\langle G \rangle_{u}$)", fontsize = 20.)

if lamb == 10.:
	plt.axis([0.0,1.0,0.0,1.4])
else:
	plt.axis([0.0,1.0,0.0,1.2])

plt.legend(loc = "upper right", frameon = False, prop ={"size":16})

plt.xticks(fontsize = 16.)
plt.yticks(fontsize = 16.)

plt.tight_layout()


script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)

if lamb == 10.:
	plt.savefig(protocell_folder + "/Figures/Getatrimorphiclambda10.png")
elif lamb == 75.:
	plt.savefig(protocell_folder + "/Figures/Getatrimorphiclambda75.png")
	





plt.show()
		
		