# Compile directly from root source in terminal.
# To compile the cython helpers modules use in terminal:
#	make all
# Then:
#	python3 Dyson_Fusion_Gas.py + params
#	sintax:
#		python3 Sim.py N Niter Nsim Nsimst beta dt lf stf
'''
	Params:
	N				->	Initial number of particles
	Niter		->	Iteration number per simulation
	Nsim		->	Number of simulations to perform
	Nsimst	->	Number label for first simulation
	beta		->	Inverse temperature
	dt			->	Virtual time step
	lf			->	fusion length
	stf			->	Iteration number to start fusion process (To include a time for thermal equilibrium)
'''
# These files need to be in the same directory

import numpy as np
import pandas as pd
import os
import argparse
import pickle
import matplotlib.pyplot as plt
from operator import itemgetter, attrgetter
import realfusionhelpers as rfh

import json
import logging
import time
import Functions as Fun

'''********************************************  Main program  ********************************************'''
inicio = time.time()

	# Takes all arguments to asign the the values of variables in the code
args = Fun.argument_parsing()

	# Set up logging info assuming loglevel is bound to the string value obtained from the command line argument
	# Convert to upper case to allow the user to specify --log=DEBUG or --log=debug
loglevel = args.log
if not args.silent:
	numeric_level = getattr(logging, loglevel.upper(), None)
	if not isinstance(numeric_level, int):
		raise ValueError('Invalid log level: %s' % loglevel)
	logging.basicConfig(format='%(asctime)s : %(message)s', level=numeric_level)

	# Name base for data files Nxx_Niterxx_Nsimxx_Nsimstxx_betaxx_dtxx_lfxx_stfxx_:
filenamebase = Fun.build_filenamebase(args)

	# Simulation parameters:
N = args.N
Niter = args.Niter
Nsim = args.Nsim
Nsimst = args.Nsimst
beta = args.beta
dt = args.dt
lf = args.lf
stf = args.stf
steps = 100

	# Folder name to save all data for parameters N, Niter, beta, lf and dt:
Folder_Name = 'N'+str(N)+'_Niter'+str(Niter)+'_Nsim'+str(Nsim)+'_beta'+str(beta)+'_dt'+str(dt)+'_lf'+str(lf)
	# Folder name to save data for each core process in the cluster (SLURM ARRAY TASK ID):
Folder_Nsimst = 'Nsimst_'+str(Nsimst)

	# Directory to save all data:
os.makedirs(Folder_Name+'/'+Folder_Nsimst, exist_ok = True)
	# Save all parameters in a json file save in Folder_Name
Fun.saveargs(args, './'+Folder_Name+'/'+Folder_Nsimst+'/Nsimst_'+str(Nsimst))

logging.info("Starting %s simulations with parameters:\n %s", Nsim,	json.dumps(vars(args), indent=4, sort_keys=True))

	# Time-evolution of frequency charges data:
Q_evol = []
	# Time evolution of N for each simulation. Allows make the statistics:
N_evol = pd.DataFrame()

	# Makes the simulation process:
for sim in range(1, Nsim+1):
		# Makes one simulation only:
	data_Ntime, data_Qtime = Fun.one_simulation_no_file(N, Niter, dt, beta, lf, stf, sim, args.freediffusion, args.Nf)

	if sim == 1:
		N_evol = data_Ntime # Since initial data frame for N_evol is empty
	else:
			# Time evolution of number particles in the system for each value of sim:
		N_evol = pd.concat([N_evol, data_Ntime], axis=1)
		# Time evolution of the system's total charge for each value of sim:
	Q_evol.append(data_Qtime)

	# Save the data on a file:
	# For N(t):
filename = './'+Folder_Name+'/'+Folder_Nsimst+'/N_t.csv'
N_evol.to_csv(filename, index = False)
	# For Q(t):
filenamehists = './'+Folder_Name+'/'+Folder_Nsimst+'/Q_t.pickle'
with open(filenamehists, 'wb') as histsfile:
	pickle.dump(Q_evol, histsfile)

logging.debug("%s", Q_evol)
logging.info("Wrote particle number results to %s", filename)
logging.info("Wrote charges histogram results to %s", filenamehists)

logging.info("Finished %s simulations", Nsim)

fin = time.time()
print('Total time for make all simulations ' + str(fin-inicio))


























