import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import numpy as np
from scipy.stats import sem
plt.rcParams["font.family"] = "Times New Roman"


lys2=[]
glu3=[]
pro4=[]
ala5=[]
pro6=[]
thr7=[]
thr8=[]
thr9=[] 
lys10=[]
glu11=[]
pro12=[]
ala13=[]
pro14=[]
thr15=[]
thr16=[]
pro17=[] 

peptide=[lys2,glu3,pro4,ala5,pro6,thr7,thr8,thr9,lys10,glu11,pro12,ala13,pro14,thr15,thr16,pro17]

lys2g=[]
glu3g=[]
pro4g=[]
ala5g=[]
pro6g=[]
thr7g=[]
thr8g=[]
thr9g=[] 
lys10g=[]
glu11g=[]
pro12g=[]
ala13g=[]
pro14g=[]
thr15g=[]
thr16g=[]
pro17g=[] 

glycopeptide=[lys2g,glu3g,pro4g,ala5g,pro6g,thr7g,thr8g,thr9g,lys10g,glu11g,pro12g,ala13g,pro14g,thr15g,thr16g,pro17g]

for res_num in range(1,17):
   res_dihedral_file = open("noglycan/rama_res"+str(res_num)+".dat")
   read_file = res_dihedral_file.readlines() 
   for i in range(0,len(read_file)):
       peptide[res_num-1].extend([read_file[i]])

for res_num in range(1,17):
   res_dihedral_file = open("withglycan/rama_res"+str(res_num)+".dat")
   read_file = res_dihedral_file.readlines()
   for i in range(0,len(read_file)):
       glycopeptide[res_num-1].extend([read_file[i]])

def structure_population(dihedrals,phi1,phi2,psi1,psi2): #list of dihedrals (in string format separated by space); ranges of phi and psi angles for a specified type of structure
   count = []
   for i in range(len(dihedrals)):
      phi, psi = dihedrals[i].split()
      if float(phi) > float(phi1) and float(phi) < float(phi2) and float(psi) > float(psi1) and float(psi) < float(psi2):
         count.extend([1])
      else: 
         count.extend([0])
   #population = np.mean(count)
   
   blocked_length = 5000 #length of block in picoseconds; will be used for calculating standard error
   blocked_averages = np.zeros(len(count)/blocked_length)

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

for i in range(len(peptide)):
  if i==5 or i==6 or i==7 or i==13 or i==14: 
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
l_alpha_population = np.array(l_alpha_population)
l_alpha_stde = np.array(l_alpha_stde)


fig, axs = plt.subplots(nrows=3,ncols=1,sharex=True) 
axs = axs.flatten()

labels = ["THR1","THR2","THR3","THR4","THR5"]
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
