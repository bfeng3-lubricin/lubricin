
dihedral_file = open("/users/bfeng3/data/bfeng3/lubricin/datasets/2repeats_peptide_glycosylated/rama_t100to200ns_300K.dat")

#Find the time points within a trajectory where every residue satisfies the PPII helical dihedral angle range
def PPII(dihedral_psi, dihedral_phi):
    count = 0
    for i in range(len(dihedral_psi)):
      if dihedral_phi[i] > -90 and dihedral_phi[i] < -20:
        if dihedral_psi[i] > 50 and dihedral_psi[i] < 240: 
          count += 1
    return count == 16

read_file = dihedral_file.readlines()
stable_time = []

#Write out the time point (in picoseconds) in the trajectory where every residue is within the PPII dihedral angle range
for i in range(0,len(read_file),16):
    current_phi = []
    current_psi = []
    for a in range(0,16):
      phi,psi = read_file[a+i].split()
      current_phi.extend([(float(phi))])
      current_psi.extend([(float(psi))])

    if(PPII(current_psi,current_phi)): 
        f = open("PPII_time.out","a")
        f.write(f"{i/16+100000}\n")
    
