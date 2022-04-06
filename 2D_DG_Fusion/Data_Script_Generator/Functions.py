import numpy as np
import pandas as pd
import os
import argparse
import json
import pickle
import time
import logging
from operator import itemgetter, attrgetter
import realfusionhelpers as rfh

def argument_parsing():
		# Parse initial parameters for the simulation:
	parser = argparse.ArgumentParser()
	parser.add_argument("N", type = int, help = "initial number of particles")
	parser.add_argument("Niter", type = int, help = "number of iterations per simulation")
	parser.add_argument("Nsim", type = int, help = "number of simulations to perform")
	parser.add_argument("Nsimst", type = int, help = "number label for first simulation")
	parser.add_argument("beta", type = float, help = "inverse temperature")
	parser.add_argument("dt", type = float, help = "virtual time step" )
	parser.add_argument("lf",type = float, help = "fusion length")
	parser.add_argument("stf", type = int, help = "iteration where to start the fusion processes")
	parser.add_argument("-s", "--silent", help = "do not print diagnose messages", action="store_true")
	parser.add_argument("--freediffusion", help = "turn off the log repulsion interaction", action="store_true")
	parser.add_argument("-Nf", type = int, help = "stops the simulation when the number of particles falls below NF. Default NF=0", default=0)
	parser.add_argument("--log", help="sets the log level. Default=INFO. DEBUG will print detailed output of the simulation", default="INFO")
	return parser.parse_args()

def build_filenamebase(args):
	'''Builds the filenamebase.
	Paramemeters:
		args: namespace with the program arguments
	Filenaming convention: 
		filename base is filename = Nxx_Niterxx_Nsimxx_Nsimstxx_betaxx_dtxx_lfxx_stfxx_
		filename + 'args.json': saves the simulations parameters
		filename + 'simXX.csv': saves the data from simulation number XX'''
	# transform arguments args to a dictionary
	args_dict = vars(args)
	# contruct the filename with the simulation parameters:
	filename = ''
	for key, val in sorted(args_dict.items()):
		if key != 'silent' and key != 'freediffusion' and key != 'Nf' and key != 'log':
			filename = filename+key+str(val)+'_'
		if args.freediffusion:
			filename = filename+'freediffusion_'
	return filename
    
def saveargs(args, filenamebase):
	# saves simulation parameters to a JSON file:
	args_dict=vars(args)
	paramsfilename = filenamebase+'_args.json'
	with open(paramsfilename, 'w') as paramsfile:
		json.dump(args_dict, paramsfile, indent=4, sort_keys=True)

def performfusions(particles, fusion_angle):
	#particles = sorted(particles, key = attrgetter('position'))
	i = 0
	while i < len(particles):
		j = i+1
		if j == len(particles):
			# i = N+1. Check fusion between particles 0 and N-1
			j = 0
		diff = np.abs(particles[i].position-particles[j].position)
		if diff > np.pi :
			diff = np.pi-diff%np.pi
		if diff < fusion_angle:
			# Fuse the two particles i and j
			particles[i].charge = particles[i].charge+particles[j].charge
			del particles[j]
		i = i+1
		# End while i
	return particles

def one_simulation(N, N_iter, dt, beta, lf, startfusion, sim, freediffusion, N_final, DT_Sim, DT_Sim_File):
	''' Performs one simulation of the Dyson brownian model on a circle with fusion events
			
			Parameters
			==========
			N : int
				Initial number of particles
			N_iter : int
				Number of iterations (steps) of the simulation
			dt : float
				Time step
			beta : float
				Inverse temperature, coupling constant
			lf : float
				Fusion length. Angle length is 0.1*lf*2*pi/N
			startfusion : int
				Iteration where to start the fusion processes
			sim : int
				Simulation number used to save the data results
			filenamebase : str
				filenamebase to save the data. The data is saved to filenamebase + srt(sim)
			N_final : int
				Final number of particles to stop the simulation if reached before N_iter '''
	
	f = 1												# Friction parameter (absorbed in unit time)
	q = 1.0											# Unit charge value
	sprom = 2*np.pi/N						# Initial system mean separation
	fusion_tol = 0.1*lf*sprom 	# Critial fusion angle
	t = 0												# Start time
	
		# To randomize first-positions of particles:
	ruido = np.random.uniform(-0.5*np.pi/N, 0.5*np.pi/N, N)
		# First positions of particles:
	firstpositions = np.linspace(0.0, 2*np.pi-2*np.pi/N, N)+ruido
	
		# Associate to each charge q its position in the circle:
	particles = [rfh.particle(q, rfh.wrapAngle(firstpositions[k])) for k in range(N)]
		# Organize all particle positions:
	particles = sorted(particles, key = attrgetter('position'))
	
		# Column names for data table of N(t):
	t_col = 't_'+str(sim)
	N_col = 'N_'+str(sim)
		# For time-evolution of Total number of particles N=N(t):
	data_n = pd.DataFrame([[t, len(particles)]], columns=[t_col, N_col])
		# Frequency data for charges:
	hist = np.histogram([each_particle.charge for each_particle in particles], np.arange(q/2, q*(N+3/2), q))
		# For time-evolution of frequency charges data:
	hists = [[sim, t, hist]]
	
	logging.info("Running simulation # %s", sim)
	logging.debug("time, N")
	logging.debug("%s, %s", t, N)
	logging.debug("%s", particles)
	logging.debug("%s", hists)
		
	k = 0
	while k < N_iter and len(particles) > N_final:
			# Computes the force and steps to move the particles:
		sigma = np.sqrt(2*dt/(f*beta))
		normalnoise = np.random.normal(0, sigma, len(particles))
		
			# Add data into DT_Sim matrix to take account of virtual time evolution of the particles positions and charges:
		DT_Sim['time'] = pd.Series([t for charges in range(len(particles))])
		DT_Sim['label'] = pd.Series([charges+1 for charges in range(len(particles))])
		DT_Sim['charge'] = pd.Series([particles[charges].charge for charges in range(len(particles))])
		DT_Sim['pos'] = pd.Series([particles[positions].position for positions in range(len(particles))])
			# Use the Langevin Eq. to find delta theta for each particle:
		rfh.stepup(particles, normalnoise, freediffusion, dt/f)
			# Add delta theta to DT_Sim matrix:
		DT_Sim['step'] = pd.Series([particles[steps].step for steps in range(len(particles))])
		
			# Move all particles by step:
		for each_particle in particles:
			each_particle.update_position()

			# Old-postition to New-position of the particles to make evolve the system:
		DT_Sim['newpos'] = pd.Series([particles[positions].position for positions in range(len(particles))])
		
			# Print data into the file: Sim_Data.dat:
		#print(DT_Sim, file = DT_Sim_File)
		if k == 0:
			DT_Sim.to_csv(DT_Sim_File, index = False)
		DT_Sim.to_csv(DT_Sim_File, header= False, index = False)
		
			# After the startfusion-th iteration check every step for fusion events
		N_now = len(particles)
		
			# Check for fusions and perform the fusions:
		if k >= startfusion:
			particles = performfusions(particles, fusion_tol)
		
			# Save only if there was a fusion (to save HDD space):
		if len(particles) < N_now:
			data_n = data_n.append({t_col: t, N_col: len(particles)}, ignore_index = True)
			hist = np.histogram([each_particle.charge for each_particle in particles], np.arange(q/2, q*(N+3/2), q))
			hists.append([sim, t, hist])
			logging.debug("%s, %s, %s", sim, t, len(particles))
			logging.debug("%s", particles)
			logging.debug("%s", hist)
		t = t+dt
		k = k+1
	# End simulation
	
	logging.info("Simulation # %s completed", sim)
	
	return data_n, hists

def one_simulation_no_file(N, N_iter, dt, beta, lf, startfusion, sim, freediffusion, N_final):
	''' Performs one simulation of the Dyson brownian model on a circle with fusion events
			Parameters
			==========
			N : int
				Initial number of particles
			N_iter : int
				Number of iterations (steps) of the simulation
			dt : float
				Time step
			beta : float
				Inverse temperature, coupling constant
			lf : float
				Fusion length. Angle length is 0.1*lf*2*pi/N
			startfusion : int
				Iteration where to start the fusion processes
			sim : int
				Simulation number used to save the data results
			filenamebase : str
				filenamebase to save the data. The data is saved to filenamebase + srt(sim)
			N_final : int
				Final number of particles to stop the simulation if reached before N_iter '''
	
	f = 1												# Friction parameter (absorbed in unit time)
	q = 1.0											# Unit charge value
	sprom = 2*np.pi/N						# Initial system mean separation
	fusion_tol = 0.1*lf*sprom 	# Critial fusion angle
	t = 0												# Start time
	
		# To randomize first-positions of particles:
	ruido = np.random.uniform(-0.5*np.pi/N, 0.5*np.pi/N, N)
		# First positions of particles:
	firstpositions = np.linspace(0.0, 2*np.pi-2*np.pi/N, N)+ruido
		# Associate to each charge q its position in the circle:
	particles = [rfh.particle(q, rfh.wrapAngle(firstpositions[k])) for k in range(N)]
		# Organize all particle positions:
	particles = sorted(particles, key = attrgetter('position'))
	
		# Column names for data table of N(t):
	t_col = 't_'+str(sim)
	N_col = 'N_'+str(sim)
		# For time-evolution of Total number of particles N=N(t):
	data_n = pd.DataFrame([[t, len(particles)]], columns=[t_col, N_col])
		# Frequency data for charges:
	hist = np.histogram([each_particle.charge for each_particle in particles], np.arange(q/2, q*(N+3/2), q))
		# For time-evolution of frequency charges data:
	hists = [[sim, t, hist]]
	
	logging.info("Running simulation # %s", sim)
	logging.debug("time, N")
	logging.debug("%s, %s", t, N)
	logging.debug("%s", particles)
	logging.debug("%s", hists)
		
	k = 0
	while k < N_iter and len(particles) > N_final:
			# Computes the force and steps to move the particles:
		sigma = np.sqrt(2*dt/(f*beta))
		normalnoise = np.random.normal(0, sigma, len(particles))

			# Use the Langevin Eq. to find delta theta for each particle:
		rfh.stepup(particles, normalnoise, freediffusion, dt/f)
		
			# Move all particles by step:
		for each_particle in particles:
			each_particle.update_position()
		
			# After the startfusion-th iteration check every step for fusion events:
		N_now = len(particles)
		
			# Check for fusions and perform the fusions:
		if k >= startfusion:
			particles = performfusions(particles, fusion_tol)
		
			# Save only if there was a fusion (to save HDD space):
		if len(particles) < N_now:
			new_data_n = pd.DataFrame([[t, len(particles)]], columns=[t_col, N_col])
			data_n = pd.concat([data_n, new_data_n], ignore_index=True)
			hist = np.histogram([each_particle.charge for each_particle in particles], np.arange(q/2, q*(N+3/2), q))
			hists.append([sim, t, hist])
			logging.debug("%s, %s, %s", sim, t, len(particles))
			logging.debug("%s", particles)
			logging.debug("%s", hist)
		t = t+dt
		k = k+1
	# End simulation
	
	logging.info("Simulation # %s completed", sim)
	
	return data_n, hists



























