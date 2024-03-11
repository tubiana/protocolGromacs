#!/bin/bash


#---------  HPC SETUP  -----------
GMX=gmx #GMX command (can be "$GMX_mpi" sometimes. Just change it here
NT=8 # Number of CPU
#THOSE COMMANDS 
GPU0="-gpu_id 0 -ntmpi 4 -ntomp 2 -nb gpu -bonded gpu -npme 1 -pme gpu -pmefft gpu -pin on -pinstride 0 -nstlist 100 -pinoffset 49" 
GPU1="-gpu_id 1 -ntmpi 4 -ntomp 2 -nb gpu -bonded gpu -npme 1 -pme gpu -pmefft gpu -pin on -pinstride 0 -nstlist 100 -pinoffset 0"
export GMX_GPU_PME_PP_COMMS=true
export GMX_FORCE_UPDATE_DEFAULT_GPU=1
export GMX_GPU_DD_COMMS=true
MDRUN_CPU="$GMX mdrun -nt ${NT}"
MDRUN_GPU="$GMX gmx mdrun $GPU0"
MDRUN="$MDRUN_CPU"
MDRUNmini=$MDRUN_CPU


source ~/miniconda3/etc/profile.d/conda.sh
conda activate gmx2023
source run_mini_only.sh
#source runGromacs.sh