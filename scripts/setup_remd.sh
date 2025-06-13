#!/bin/bash
#This script sets up each replica directory. 

i=1
module load gromacs/2018.2
while read -r line; do
  mkdir rep$i
  cd rep$i 
  cp ../* .
  cp -r ../toppar .
  sed -i 's/TEMP/'"$line"'/g' nvt.mdp
  sed -i 's/TEMP/'"$line"'/g' npt.mdp
  sed -i 's/TEMP/'"$line"'/g' md.mdp
  gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
  i=$((i + 1))
  cd ..
done < temperature.dat
