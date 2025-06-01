#!/bin/bash
#SBATCH -n 128
#SBATCH --mem-per-cpu=4G
#SBATCH -t 1:00:00
#SBATCH -J gromacs
#SBATCH --mail-user=bibo_feng@brown.edu


module load gromacs/2018.2_hpcx_2.7.0_gcc_10.2_slurm20
module load mpi/hpcx_2.7.0_gcc_10.2_slurm20 gcc/10.2

 
srun --mpi=pmix gmx_mpi mdrun -deffnm nvt -multidir rep1 rep2 rep3 rep4 rep5 rep6 rep7 rep8 rep9 rep10 rep11 rep12 rep13 rep14 rep15 rep16 rep17 rep18 rep19 rep20 rep21 rep22 rep23 rep24 rep25 rep26 rep27 rep28 rep29 rep30 rep31 rep32 rep33 rep34 rep35 rep36 rep37 rep38 rep39 rep40 rep41 rep42 rep43 rep44 rep45 rep46 rep47 rep48 rep49 rep50 rep51 rep52 rep53 rep54 rep55 rep56 rep57 rep58 rep59 rep60 rep61 rep62 rep63 rep64 -maxh 48

