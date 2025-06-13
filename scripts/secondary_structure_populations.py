import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import numpy as np
from scipy.stats import sem
plt.rcParams["font.family"] = 'serif' 

only_plot_threonines = True #Switch to false if you want to plot secondary structure populations of every residue

glycopeptide_rama_directory = "/users/bfeng3/data/bfeng3/lubricin/datasets/2repeats_peptide_glycosylated/" 
peptide_rama_directory = "/users/bfeng3/data/bfeng3/lubricin/datasets/2repeats_peptide_noglycan/" 

num_sampled_dihedrals = 100000 #dihedral angles are sampled and recorded every 100 ns / 100000 = 1 ps  
dihedral_labels = ["LYS2","GLU3","PRO4","ALA5","PRO6","THR7","THR8","THR9","LYS10","GLU11","PRO12","ALA13","PRO14","THR15","THR16","PRO17"]

#2D array to store pairs of dihedral angles for each residue at every sampled time point
peptide = np.empty((len(dihedral_labels),num_sampled_dihedrals),dtype=object)
glycopeptide = np.empty((len(dihedral_labels),num_sampled_dihedrals),dtype=object)


for res_num in range(len(dihedral_labels)):
   peptide_dihedral_file = open(f"{peptide_rama_directory}/rama_"+f"{dihedral_labels[res_num]}"+"_t100to200ns_300K.dat")
   peptide_read_file = peptide_dihedral_file.readlines() 
   
   glycopeptide_dihedral_file = open(f"{glycopeptide_rama_directory}/rama_"+f"{dihedral_labels[res_num]}"+"_t100to200ns_300K.dat")
   glycopeptide_read_file = glycopeptide_dihedral_file.readlines()
  
   for i in range(0,num_sampled_dihedrals):
       #peptide[res_num][i] or glycopeptide[res_num][i] has the format ["phi psi"]
       peptide[res_num][i] = peptide_read_file[i] 
       glycopeptide[res_num][i] = glycopeptide_read_file[i]


#dihedrals = array of dihedral angles (in string format separated by space); ranges of phi and psi angles for a specified type of structure
def structure_population(dihedrals,phi1,phi2,psi1,psi2): 
   count = []
   for i in range(len(dihedrals)):
      phi, psi = dihedrals[i].split()
      if float(phi) > float(phi1) and float(phi) < float(phi2) and float(psi) > float(psi1) and float(psi) < float(psi2):
         count.extend([1])
      else: 
         count.extend([0])
   
   blocked_length = 5000 #length of block in picoseconds; will be used for calculating standard error
   blocked_averages = np.zeros(int(len(count)/blocked_length))

   for i in range(len(blocked_averages)):
       blocked_averages[i] = sum(count[i*blocked_length:(i+1)*blocked_length])/float(blocked_length)
   population = np.mean(blocked_averages)    
   stderror = 1.96*sem(blocked_averages) 
   return [population,stderror] 

PPII_population = []
PPII_stde = [] #PPII standard error
beta_population = []
beta_stde = []
r_alpha_population = []
r_alpha_stde = []
l_alpha_population = []
l_alpha_stde = []


if only_plot_threonines: 
    residues_to_plot = [5,6,7,13,14] #These are the residue indices for threonines 
else:
    residues_to_plot = [i for i in range(len(dihedral_labels))]
 
for i in residues_to_plot:
     PPII_population.append([structure_population(peptide[i],-90,-20,50,240)[0],structure_population(glycopeptide[i],-90,-20,50,240)[0]])
     PPII_stde.append([structure_population(peptide[i],-90,-20,50,240)[1],structure_population(glycopeptide[i],-90,-20,50,240)[1]])
     beta_population.append([structure_population(peptide[i],-180,-90,50,240)[0],structure_population(glycopeptide[i],-180,-90,50,240)[0]])
     beta_stde.append([structure_population(peptide[i],-180,-90,50,240)[1],structure_population(glycopeptide[i],-180,-90,50,240)[1]])
     r_alpha_population.append([structure_population(peptide[i],-160,-20,-120,50)[0],structure_population(glycopeptide[i],-160,-20,-120,50)[0]])
     r_alpha_stde.append([structure_population(peptide[i],-160,-20,-120,50)[1],structure_population(glycopeptide[i],-160,-20,-120,50)[1]])


PPII_population = np.array(PPII_population)
PPII_stde = np.array(PPII_stde)
beta_population = np.array(beta_population)
beta_stde = np.array(beta_stde)
r_alpha_population = np.array(r_alpha_population)
r_alpha_stde = np.array(r_alpha_stde)

#Plot Secondary Structure Populations 
fig, axs = plt.subplots(nrows=3,ncols=1,sharex=True) 
axs = axs.flatten()

if only_plot_threonines:
    labels = ["THR1","THR2","THR3","THR4","THR5"]
else:
    labels = dihedral_labels

x = np.arange(1,len(labels)+1,1)

width=0.125

axs[0].bar(x-width/2,PPII_population[:,0],width,label="Peptide",yerr=PPII_stde[:,0],color="purple",capsize=3)
axs[0].bar(x+width/2,PPII_population[:,1],width,label="Glycopeptide",yerr=PPII_stde[:,1],color="red",capsize=3)
axs[1].bar(x-width/2,beta_population[:,0],width,label="Peptide",yerr=beta_stde[:,0],color="purple",capsize=3)
axs[1].bar(x+width/2,beta_population[:,1],width,label="Glycopeptide",yerr=beta_stde[:,1],color="red",capsize=3)
axs[2].bar(x-width/2,r_alpha_population[:,0],width,yerr=r_alpha_stde[:,0],color="purple",capsize=3)
axs[2].bar(x+width/2,r_alpha_population[:,1],width,yerr=r_alpha_stde[:,1],color="red",capsize=3)

plt.xticks(x,labels,rotation='vertical')
axs[0].set_ylabel("PPII")#,fontweight="bold")
axs[1].set_ylabel(r'$\beta$' ,fontweight="bold")
axs[2].set_ylabel(r'$\alpha_R$' ,fontweight="bold")

axs[0].set_ylim(0,1)
axs[1].set_ylim(0,1)
axs[2].set_ylim(0,1)

plt.xlabel("Residue",fontweight="bold")
axs[1].legend()
plt.savefig("PPII_population_change.png",dpi=600)
plt.show()
