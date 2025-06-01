#!/bin/bash

i=1
module load gromacs/2018.2
while read -r line; do
  cd rep$i 
  gmx grompp -f md.mdp -c npt.gro -t npt.cpt  -p topol.top -n index.ndx -o md.tpr
  i=$((i + 1))
  cd ..
done < temperature.dat
