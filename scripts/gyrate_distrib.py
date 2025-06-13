#This script reads output file of gmx gyrate to calculate radius of gyration distribution across a trajectory, after discarding the initial {equilibration_time} ps of the trajectory 
import matplotlib.pyplot as plt
import numpy as np 


equilibration_time = 100000 #in picoseconds 
num_bins = 50
filename = "../datasets/2repeats_peptide_glycosylated/gyrate_300K.xvg" #gmx gyrate output file. Delete the lines before the actual data in the initial output file before using this script
plot_title = "Glycopeptide"

#opens and reads in {filename} file , which should contain Rg (nm) at every picosecond of trajectory 
Rg_file = open(filename,"r+") 
Rg = [] 

#loops through lines of Rg_file to add radius to Rg list
for line in Rg_file: 
  #split each line by space; in a typical GROMACS output file, the second column is what contains the Rg of interest
  line_split = line.split()
  Rg.append(float(line_split[1]))

#Disgards the first {equilibration_time} picoseconds of trajectory by removing the first {equilibration_time} items in Rg list 
for i in range(0,equilibration_time):
  Rg.pop(0) 

#Find midpoint of the remaining trajectory, then split Rg array into two halves   
midpoint = int(len(Rg)/2)
Rg1 = Rg[0:midpoint]  #list of Rg in 1st time interval (between {equilibration_time} picoseconds and {midpoint} picoseconds)
Rg2 = Rg[midpoint:]  #list of Rg in 2nd time interval (after {midpoint} picoseconds)


#calculate the size of each histogram bin, which is just the range of radius divided by total number of bins
minimum_Rg = min(Rg)
maximum_Rg = max(Rg)
Rg_range = maximum_Rg - minimum_Rg
bin_size = Rg_range/num_bins 

#create list of histogram bins
bins = [] 
i = 0
while i < num_bins:
  bins.append(minimum_Rg + i*bin_size) 
  i += 1

#construct probability_Rg1 and probability_Rg2, which are arrays that contain the probability of falling into each histogram bin during the first or second time interval
probability_Rg1 = np.zeros(num_bins)  
probability_Rg2 = np.zeros(num_bins)
for bin_index in range(num_bins): 
  #Loop through Rg1 and Rg2 to count the number of times a radius falls into each bin during each of the two time intervals
  for i1 in range(len(Rg1)): 
    if Rg1[i1] >= bins[bin_index] and Rg1[i1] < bins[bin_index] + bin_size:
      probability_Rg1[bin_index] += 1

  for i2 in range(len(Rg2)): 
    if Rg2[i2] >= bins[bin_index] and Rg2[i2] < bins[bin_index] + bin_size:
      probability_Rg2[bin_index] += 1

#convert to probabilities by dividing list of counts by total number of data points at each time interval   
probability_Rg1 = probability_Rg1/len(Rg1)
probability_Rg2 = probability_Rg2/len(Rg2) 


#calculate the proportion of overlap of each distribution by summing up the distribution with the lower probability at each bin
overlap_proportion = 0
for i in range(num_bins):
  overlap_proportion += min(probability_Rg1[i], probability_Rg2[i])  

if overlap_proportion > .90: 
  print(f"The percentage of overlap is {round(overlap_proportion*100,2)}%, and the system is well converged!")
else:
  print(f"The percentage of overlap is only {round(overlap_proportion*100,2)}%, so the system has not converged.")
  

#create a figure of the approximated probability distributions by plotting "bins" against "probability_Rg[1 or 2]"  
fig, ax = plt.subplots()
label1 = "{} ns to {} ns".format(int(equilibration_time/1000), int((midpoint+equilibration_time)/1000)) #label for the starting and ending timepoints of the first time interval (divided by 1000 to convert to nanoseconds)
label2 = "{} ns to {} ns".format(int((midpoint+equilibration_time)/1000), int((len(Rg)+equilibration_time)/1000)) #label for the starting and ending timepoints of the second time interval (divided by 1000 to convert to nanoseconds)
ax.plot(bins,probability_Rg1,label=label1)
ax.plot(bins,probability_Rg2,label=label2)
ax.set_ylabel("Density")
ax.set_xlabel("Radius of Gyration (nm)")
ax.legend()
ax.text(0.8,0.04,"Percentage Overlap: {}%".format(round(overlap_proportion*100,2)))
plt.title(f"{plot_title}")
plt.savefig("gyrate_distribution.png",dpi=600)
plt.show()


