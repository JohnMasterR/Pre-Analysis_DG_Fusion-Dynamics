# Compile directly from root source in terminal.
#	python3 Q_t_Analysis.py

import numpy as np
import pandas as pd
import pickle

#+++++++++++++++++Functions+++++++++++++++++#

def time_df(dt, N_iter, steps):
    '''Generates time dataframe (x-axis) with steps uniformly on a log scale
    Parameters:
        dt: float    -> Original time step
        N_iter: int  -> Total number of iterations (time steps) in the simulation
        steps: int   -> number of time steps to generate 
    Returns:
        t_df: pd.Dataframe	-> The time dataframe'''
    tf = N_iter*dt
    t_df = pd.DataFrame([(t) for t in np.geomspace(dt, tf, steps)], columns = ['t'])
    return t_df

def empty_Q_hist(Niter, Nsim, dt, N, Qmax, steps):
    Charge_Cols = [] # Column names for stat. analysis of charge distributions
    for q in range(Qmax):
        Charge_Cols.append('Q'+str(q+1))
    hists_out = time_df(dt, Niter, steps)
    Cols_Charges = pd.DataFrame(np.zeros((len(hists_out), Qmax)), columns=Charge_Cols)

    hist_out = pd.concat([hists_out, Cols_Charges], axis=1)
    return hist_out

def empty_Q_pd(Niter, Nsim, dt, N, Qmax, steps):
    Charge_Sim = [] # Column names for stat. analysis of charge distributions
    for sim in range(Nsim):
        for q in range(Qmax):
            Charge_Sim.append('Sim'+str(sim+1)+'Q'+str(q+1))

    hists_out = time_df(dt, Niter, steps)
    Cols_Charges = pd.DataFrame(np.zeros((len(hists_out), Nsim*Qmax)), columns=Charge_Sim)

    hist_out = pd.concat([hists_out, Cols_Charges], axis=1)
    return hist_out

def decompress_data_Q(compressed_hists, N, Qmax, Niter, Nsim, dt, steps):
    hist_out = empty_Q_pd(Niter, Nsim, dt, N, Qmax, steps)
    for sim in range(Nsim):
        for i, line in enumerate(compressed_hists[sim]):
            s, t_now, (hist, bins) = compressed_hists[sim][i]
            if i < len(compressed_hists[sim])-1:
                tf = compressed_hists[sim][i+1][1]
            else:
                tf = Niter*dt
            for idx, Q in enumerate(hist):
                col_name = 'Sim'+str(sim+1)+'Q'+str(idx+1)
                hist_out.loc[(t_now <= hist_out['t']) & (hist_out['t'] <= tf), [col_name]] = Q
    return hist_out

#+++++++++++++++++Parameters+++++++++++++++++#

N = 100
Nsim = 10
beta = 0.5
lf = 1.0

steps = 100
Qmax = N

Arr_Used = np.arange(0,100) # Number of arrays used in sbatch: 0 - 99

#+++++++++++++++++Read and extract interest data for Q(t)+++++++++++++++++#

dt1 = 1e-4
Niter1 = 1000000
directory1 = 'N'+str(N)+'_Niter'+str(Niter1)+'_Nsim'+str(Nsim)+'_beta'+str(beta)+'_dt'+str(dt1)+'_lf'+str(lf)

DT = pd.DataFrame()
for i in range(5):#Arr_Used:
	print('Q_t start for sim = '+str(i+1))
	hist_filename = directory1+'/Nsimst_'+str(i)+'/Q_t.pickle'
	with open(hist_filename, 'rb') as histsfile:
		hists=pickle.load(histsfile)

	tmp_df = decompress_data_Q(hists, N, Qmax, Niter1, Nsim, dt1, steps)
	DT = pd.concat([DT, tmp_df], axis=1)
	print('Q_t end for sim = '+str(i+1))

DT.to_csv(directory1+'_All_Q.csv')



#Q_hist_ave = empty_Q_hist(Niter, Nsim, dt, N, Qmax, steps)
#Q_hist_std = empty_Q_hist(Niter, Nsim, dt, N, Qmax, steps)


























