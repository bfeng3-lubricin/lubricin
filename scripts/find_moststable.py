def PPII(dihedral_psi, dihedral_phi):
    count = 0
    for i in range(len(dihedral_psi)):
      if dihedral_phi[i] > -90 and dihedral_phi[i] < 20:
        if dihedral_psi[i] > 50 and dihedral_psi[i] < 240: 
          count += 1
    return count == 16
dihedral_file = open("rama.dat")
read_file = dihedral_file.readlines()
stable_time = []
for i in range(0,3200016,16):
    current_phi = []
    current_psi = []
    for a in range(0,16):
      phi,psi = read_file[a+i].split()
      current_phi.extend([(float(phi))])
      current_psi.extend([(float(psi))])
    if(PPII(current_phi,current_psi)): 
        stable_time.extend([i/16])
print(stable_time)
