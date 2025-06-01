#!/bin/bash

i=1
module load gromacs/2018.2
while read -r line; do
  cd rep$i 
  gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
  i=$((i + 1))
  cd ..
done < temperature.dat
