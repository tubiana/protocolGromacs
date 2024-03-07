#!/bin/bash
#SBATCH --partition=cpu_p1
#SBATCH --job-name="TT:MD"
#SBATCH --nodes=8
#SBATCH --ntasks=320
#SBATCH --ntasks-per-node=40
#SBATCH --hint=nomultithread # pas d'hyperthreading
#SBATCH --time=19:55:00 # Temps dexecution max
#SBATCH --output=Hybride%j.out # fichier de sortie
#SBATCH --error=Hybride%j.out # fichier d'erreur
#SBATCH --account=zvk@cpu


# on se place dans le répertoire de soumission
cd ${SLURM_SUBMIT_DIR}

module purge
module load cpuarch/amd
module load gromacs/2023.4-mpi-cuda

#Author: Thibault Tubiana, PhD, 2020
#Version: 1.3.0




#---------  HPC SETUP  -----------
MPI="" #If you have to submit jobs with MPI softwares like "mpirun -np 10". Add the command here
GMX=gmx_mpi #GMX command (can be "$GMX_mpi" sometimes. Just change it here
#THOSE COMMANDS 
GPU0="-gpu_id 0 -ntmpi 4 -ntomp 2 -nb gpu -bonded gpu -npme 1 -pme gpu -pmefft gpu -pin on -pinstride 0 -nstlist 100 -pinoffset 49" 
GPU1="-gpu_id 1 -ntmpi 4 -ntomp 2 -nb gpu -bonded gpu -npme 1 -pme gpu -pmefft gpu -pin on -pinstride 0 -nstlist 100 -pinoffset 0"
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=1
export GMX_GPU_DD_COMMS=true
MDRUN_CPU="srun $GMX mdrun"
MDRUN_GPU="$GMX gmx mdrun -ntomp 2 -ntmpi 4 -npme 1 -bonded gpu -nb gpu -pme gpu -pmefft gpu"
MDRUN="$MDRUN_CPU"
MDRUNmini=$MDRUN_CPU

bash runGromacs.sh