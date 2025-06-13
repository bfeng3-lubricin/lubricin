#This script generates a Ramachandran free energy surface given a file containing a list of dihedral angles sampled in an MD simulation
#Input: file containing list of dihedral angles sampled in MD subjectory 

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

filename = sys.argv[1] #enter file containing list of dihedral angles sampled in MD trajectory. Example input file can be found in ../datasets/2repeats_peptide_noglycan/rama_t100to200ns_300K.dat 

#Relevant Physical Constants and Parameters 
T = 300 #Temperature in Kelvin
kB = 3.2976268E-27 #Boltzmann Constant in kcal/K
Avogadro = 6.02214179E23 #Avogadro's Number

#Creating necessary arrays
phi = np.linspace(-180,180,91) #results in bins of 4 degrees x 4 degrees
psi = np.linspace(-180,180,91)
Psi_grid, Phi_grid = np.meshgrid(phi,psi) 
grid_interval = phi[-1] - phi[-2] #size of histogram bins (4 degrees in this case)
count = np.zeros((91,91))
dG = np.zeros((91,91))

#Reads file containing list of dihedral angles sampled in MD trajectory
open_rama = open(filename,"r")
read_rama = open_rama.readlines()

#Binning Dihedral Angles
for i in range(len(read_rama)):
  phi_data = float(read_rama[i].split()[0]) 
  phi_data = phi_data - phi_data % grid_interval #Rounds down to nearest bin 
  phi_index = np.where(phi == phi_data)[0]

  psi_data = float(read_rama[i].split()[1])
  psi_data = psi_data - psi_data % grid_interval
  psi_index = np.where(psi == psi_data)[0]

  count[phi_index[0]][psi_index[0]] += 1


#Calculating Free Energy for Each Bin 
max_count = max(count.flatten()) #Will be used for normalization, such that the lowest free energy = 0

for x in range(len(phi)):
  for y in range(len(psi)):
    if count[x][y] != 0:
        dG[x][y] = -1*kB*T*Avogadro*np.log(count[x][y]/max_count) #free energy in kcal/mol 
    else:
        dG[x][y] = float("inf")

#Plots Free Energy Surface
plt.contourf(Phi_grid,Psi_grid,dG,cmap=cm.gnuplot2,vmax=6.4)
plt.xlim(-180.180) 
plt.ylim(-180,180)
plt.xticks(ticks=[-180,-120,-60,0,60,120,180],labels=[-180,-120,-60,0,60,120,180])
plt.yticks(ticks=[-180,-120,-60,0,60,120,180],labels=[-180,-120,-60,0,60,120,180])
plt.xlabel(r'$\phi$ (Degrees)')
plt.ylabel(r'$\psi$ (Degrees)')
cbar = plt.colorbar()
cbar.set_label(r'$\Delta G$'+' (kcal/mol)',size=12)
plt.savefig(filename.split(".")[0]+str("_FES.png"),dpi=600)
plt.show()
