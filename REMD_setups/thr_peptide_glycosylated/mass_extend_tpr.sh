#!/bin/bash

i=1
module load gromacs/2018.2
while read -r line; do
  cd rep$i 
  gmx convert-tpr -s extend_md.tpr -extend 1300000 -o md.tpr 
  i=$((i + 1))
  cd ..
done < temperature.dat
