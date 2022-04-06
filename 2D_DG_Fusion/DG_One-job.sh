#!/bin/bash

# usage: DG_One-job.py N Niter Nsim Nsimst beta dt lf stf [-h] [-s] [--freediffusion] [-Nf NF] [--log LOG]
#

#SBATCH --job-name=Test-fusion-script	#Nombre del job
#SBATCH -p medium			#Cola a usar, Default=short (Ver colas y límites en /hpcfs/shared/README/partitions.txt)
#SBATCH -N 1				#Nodos requeridos, Default=1
#SBATCH -n 1				#Tasks paralelos, recomendado para MPI, Default=1
#SBATCH --cpus-per-task=1		#Cores requeridos por task, recomendado para multi-thread, Default=1
#SBATCH --mem=2048			#Memoria en Mb por CPU, Default=2048
#SBATCH --time=7-0:00:00		#Tiempo maximo de corrida, Default=2 horas
#SBATCH -o Test-fusion-script.o%A_%a	#Nombre de archivo de salida

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
Niter=10000000
Nsim=1
Nsimst=1
stf=0

module load anaconda/python3.9
conda init bash
source ~/bin/conda_init.sh
conda activate cythonenv

for lf in 0.1 1.0 4.0 7.0
do
	for beta in 0.1 1.0 3.0
	do
		for dt in 1e-6 1e-7 1e-8
		do
			inicio=`date +%s%N`
			python3 Data_Script_Generator/Dyson_Fusion_Gas.py $N $Niter $Nsim $Nsimst $beta $dt $lf $stf &
			fin=`date +%s%N`
			T_time=$(($fin-$inicio))
			printf "$lf, $beta, $dt, $T_time\n" >> Times_Niter10000000_CPU99.csv
		done
	done
done
wait

