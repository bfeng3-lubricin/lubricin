#!/bin/bash
#SBATCH --ntasks-per-node 32
#SBATCH -N 1
#SBATCH --mem-per-cpu=4G
#SBATCH -t 1:00:00
#SBATCH -J gromacs
#SBATCH --mail-user=bibo_feng@brown.edu


module load gromacs/2018.2_hpcx_2.7.0_gcc_10.2_slurm20
module load mpi/hpcx_2.7.0_gcc_10.2_slurm20 gcc/10.2
 
srun --mpi=pmix gmx_mpi mdrun -deffnm nvt -multidir rep1 rep2 rep3 rep4 rep5 rep6 rep7 rep8 rep9 rep10 rep11 rep12 rep13 rep14 rep15 rep16 rep17 rep18 rep19 rep20 rep21 rep22 rep23 rep24 rep25 rep26 rep27 rep28 rep29 rep30 rep31 rep32 -maxh 1

