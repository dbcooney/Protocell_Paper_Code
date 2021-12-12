"""
Script used to generate Figure 3.1, illustrating the phase portrait and sample
trajectories for the within-protocell competition between fast, slow, and dimer 
replicators.
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
bS = 1.
bF = 1. + s
bD = (1. + s) / (2. + s)

timesteps = 8000
time_diff = 0.01

x_step = 0.001

N = 20

"""
Righthand sides of ODEs for within-protocell competition.
"""

def xrighthand(x,y,bS,bF,bD):
	return x * (bS - bD + (bD - bS) * x + (bD - bF) * y)
	
def yrighthand(x,y,bS,bF,bD):
	return y * (bF - bD + (bD - bS) * x + (bD - bF) * y)
	
def simplex_edge(x):
	return 1.0 - x
	
state1_init = [0.1,0.1]
state1 = state1_init 

state2_init = [0.31,0.01]
state3_init = [0.51,0.004]

statex = [state1_init[0]]
statey = [state1_init[1]]

statex2 = [state2_init[0]]
statey2 = [state2_init[1]]

statex3 = [state3_init[0]]
statey3 = [state3_init[1]]




"""
Integrating sample within-protocell trajectories in time.
"""

for time in range(timesteps):
	
	
	xold = statex[-1]
	yold = statey[-1]
	
	xnew = xold + time_diff * xrighthand(xold,yold,bS,bF,bD)
	ynew = yold + time_diff * yrighthand(xold,yold,bS,bF,bD)
	
	statex.append(xnew)
	statey.append(ynew)
	
	
	x2old = statex2[-1]
	y2old = statey2[-1]
	
	x2new = x2old + time_diff * xrighthand(x2old,y2old,bS,bF,bD)
	y2new = y2old + time_diff * yrighthand(x2old,y2old,bS,bF,bD)
	
	statex2.append(x2new)
	statey2.append(y2new)
	
	x3old = statex3[-1]
	y3old = statey3[-1]
	
	x3new = x3old + time_diff * xrighthand(x3old,y3old,bS,bF,bD)
	y3new = y3old + time_diff * yrighthand(x3old,y3old,bS,bF,bD)
	
	statex3.append(x3new)
	statey3.append(y3new)
	


x_vec = np.arange(0.,1. + x_step,x_step)



"""
Calculating and plotting arrows for within-protocell phase portrait.
"""
for i in range(N):
	for j in range(N):
		x = np.float(i)/N
		y = np.float(j)/N
		
		dx = xrighthand(x,y,bS,bF,bD)
		dy = yrighthand(x,y,bS,bF,bD)
		
		if i + j <= N:
			plt.quiver(x,y,dx,dy, color = 'b', width = 0.005, alpha = 0.8)
	
"""
Plotting sample within-protocell trajectories.
"""		
plt.plot(statex,statey, lw = 6., color = 'r', alpha = 0.7)	
plt.plot(statex2,statey2, lw = 6., color = 'r', alpha = 0.7, ls = '--')	
plt.plot(statex3,statey3, lw = 6., color = 'r', alpha = 0.7, ls = '-.')	
plt.plot(x_vec,simplex_edge(x_vec), color = 'k', lw = 4., alpha = 0.7)		

plt.tick_params(top = False, right = False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)

plt.axis([0.0,1.,0.0,1.])

plt.xticks(fontsize = 14.)
plt.yticks(fontsize = 14.)

plt.xlabel(r"Fraction of Slow Replicators $x$", fontsize = 20.,labelpad = 10.)
plt.ylabel(r"Fraction of Fast Replicators $y$", fontsize = 20.)

script_folder = os.getcwd()
protocell_folder = os.path.dirname(script_folder)
plt.savefig(protocell_folder + "/Figures/withinprotocellphase.png",transparent = True,bbox_inches='tight',pad = 0)


plt.show()
	
	