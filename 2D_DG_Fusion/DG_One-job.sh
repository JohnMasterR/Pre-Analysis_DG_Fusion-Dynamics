#!/bin/bash

# usage: Data_Script_Generator/Dyson_Fusion_Gas.py N Niter Nsim Nsimst beta dt lf stf
#                        				[-h] [-s] [--freediffusion] [-Nf NF] [--log LOG]

#SBATCH --array=0-99	#Lanza 100 Jobs numerados de 0 hasta 100
#SBATCH --job-name=real_fusion_job		#Nombre del job
#SBATCH -p medium		#Cola a usar, Default=short (Ver colas y límites en /hpcfs/shared/README/partitions.txt)
#SBATCH -N 1				#Nodos requeridos, Default=1
#SBATCH -n 1				#Tasks paralelos, recomendado para MPI, Default=1
#SBATCH --cpus-per-task=1		#Cores requeridos por task, recomendado para multi-thread, Default=1
#SBATCH --mem=2048		#Memoria en Mb por CPU, Default=2048
#SBATCH --time=7-0:00:00			#Tiempo máximo de corrida, Default=2 horas
#SBATCH -o N100_Niter100000000_Nsim10_beta0.5_dt1e-06_lf1.0/real_fusion_job.o%A_%a		#Nombre de archivo de salida

# positional arguments:
#   N                initial number of particles
#   Niter            number of iterations per simulation
#   Nsim             number of simulations to perform
#   Nsimst           number label for first simulation
#   beta             inverse temperature
#   dt               virtual time step
#   lf               fusion length
#   stf              iteration where to start the fusion processes

# optional arguments:
#   -h, --help       show this help message and exit
#   -s, --silent     do not print diagnose messages
#   --freediffusion  turn off the log repulsion interaction
#   -Nf NF           stops the simulation when the number of particles falls
#                    below NF. Default NF=0
#   --log LOG        sets the log level. Default=INFO. DEBUG will print
#                    detailled output of the simulation

N=100
Niter=100000000
Nsim=10
Nsimst=$(( SLURM_ARRAY_TASK_ID ))
beta=0.5
dt=1e-6
lf=1.0
stf=0

module load anaconda/python3.9

python3 Data_Script_Generator/Dyson_Fusion_Gas.py $N $Niter $Nsim $Nsimst $beta $dt $lf $stf 

