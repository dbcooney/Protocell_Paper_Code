Code for "A PDE Model for Protocell Evolution and the Origin of Chromosomes via Multilevel Selection"
==============

This repository accompanies the preprint "A PDE Model for Protocell Evolution and 
the Origin of Chromosomes via Multilevel Selection", by Daniel B. Cooney, 
Fernando W. Rossine, Dylan H. Morris, and Simon A. Levin. It includes all of the scripts
used to compute numerical finite volume simulations of the dimorphic and trimorphic
multilevel selection models and to generate all of the figures in the paper. The
repository also includes all of the figures from the paper, as well as the data 
produced by numerical simulations that require the use of additional scripts for 
generating figures. 

The repository is organized into the following folders: Scripts, Simulation_Outputs, and
Figures. All of the scripts can be run using Python 2.7.

The main Python scripts used to compute numerical finite volume simulations are 
fv_protocell.py (for the dimorphic multilevel competition) and fvfsdtrimorphicdynamics.py
(for the trimorphic dynamics). Additional simulations were run in the dimorphic case to 
compare with analytically calculated steady states in etasteadyfv.py, and additional
simulations of the trimorphic multilevel dynamics were run in loopfvtrimorphic.py
and etaloopfvtrimorphic.py to generate data for numerical solutions as a function of the
relative strength of between-protocell selection and the complementarity of the fast
and slow genes, respectively.

For reference, below is a list of figures and the scripts and simulation output files 
that were used to generate each figure.

- Figure 2.1: Run Gofxfs.py
- Figure 2.2: Run etamodeldensities.py
- Figure 2.3: Run etapeak.py
- Figure 3.1: Run withintrimorphic.py
- Figure 3.2: Run fvfsdtrimorphicdynamics.py for each value of the parameter eta
- Figure 4.1: Run Gofzfd.py
- Figure 4.2: Run fastdimerdensities.py
- Figure 4.3: Run fvprotocell.py with FS and FD options
- Figure 4.4: Run FSvsFDprotocellfitness.py
- Figure 4.5: Run FSvsFDlargelambdacomparison.py
- Figure 4.6: Run etasregimes.py
- Figure 4.7: Run Gofzsd.py
- Figure 4.8: Run slowdimerpeakplot.py
- Figure 5.1: Run fvfsdtrimorphicdynamics.py with "quantity = trajectories"
- Figure 5.2: Run fvfsdtrimorphicdynamics.py to generate data, then run avgGintime.py
- Figure 5.3: Run fvfsdtrimorphicdynamics.py with "quantity = steady"
- Figure 5.4: Run loopfvtrimorphic.py to generate data, then run Gplottrimorphic.py
- Figure 5.5: Run loopfvtrimorphic.py to generate data, then run peakmeantrimorphic.py
- Figure 5.6: Run etaloopfvtrimorphic.py
- Figure B.1: Run etasteadyfv.py with "edge = FS"
- Figure B.2: Run etasteadyfv.py with "edge = FD"
			  