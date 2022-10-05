from biopandas.pdb import PandasPdb
import pandas as pd
import os
import re

# Antechamber Ligand Treatment
os.system("antechamber -i ligand.mol2 -fi mol2 -o ligand.in -fo prepi -c bcc -nc 0")
os.system("parmchk2 -i ligand.in -o ligand.frcmod -f prepi -a Y")

# Deleting unneeded files
os.system("rm -rf ANTECHAMBER*")
os.system("rm -rf leap.log")
os.system("rm -rf sqm*")
os.system("rm -rf ATOMTYPE.INF PREP.INF")
os.system("mv NEWPDB.PDB ligand_antechambered.pdb")

salt_conc = 150
resids = ["ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLN", "GLU", "GLX", "GLY", "HIS", "ILE", 
          "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ACE", "Y2P", 
          "HIE", "NME", "HID", "F9J", "WAT"]
weights_Da = [89, 174, 132, 133, 133, 121, 146, 147, 147, 75, 155, 131, 131, 146, 149, 165,
              115, 105, 119, 204, 181, 117, 43, 261, 155, 31, 155, 0, 0]
resids_weights = pd.DataFrame({'resids': resids, 'weights_Da': weights_Da})

with open('input_salt.leap', 'w') as f:
    f.write('''
loadamberparams Y2P.frcmod
loadOff Y2P.off
loadamberprep ligand.in
loadamberparams ligand.frcmod
source leaprc.protein.ff14SB
source leaprc.gaff
source leaprc.water.tip4pew
set default FlexibleWater on
set default PBRadii mbondi2
WAT = T4E
HOH = T4E
loadAmberParams frcmod.ionsjc_tip4pew
loadAmberParams frcmod.tip4pew
pdb = loadpdb system.pdb
charge pdb
solvateOct pdb TIP4PEWBOX 10
savepdb pdb system_salt.pdb
quit
''')
os.system("tleap -f input_salt.leap > input_salt.out")
ppdb = PandasPdb()
ppdb.read_pdb("system_salt.pdb")
df = ppdb.df["ATOM"]
df_grouped = df.groupby('residue_name')['residue_number'].nunique()
pdb_resid_list = df_grouped.index.to_list()
#print(pdb_resid_list)
pdb_resid_frequency_list = df_grouped.values.tolist()
#print(pdb_resid_frequency_list)
weights_list = []
for i in pdb_resid_list:
    weights_list.append(resids_weights["weights_Da"][resids_weights.loc[resids_weights.isin([i]).any(axis=1)].index.tolist()[0]])
#print(weights_list)
protein_weight = sum([a*b for a,b in zip(weights_list,pdb_resid_frequency_list)])/1000
#print(protein_weight)
with open("system_salt.pdb","r") as f:
    water_lines = []   
    for line in f:
        if "EPW" in line:
            water_lines.append(line)
no_water = len(water_lines)
with open("input_salt.out","r") as file: 
    lines = file.readlines()
    for line in lines:
        if line.startswith("Total perturbed charge:"):
            charge = float(round(float(re.findall(r'[-+]?\d+[,.]?\d*',line)[0])))
os.system("rm -rf input_salt.leap system_salt.pdb leap.log input_salt.out")
print("Please go the website " , "https://www.phys.ksu.edu/personal/schmit/SLTCAP/SLTCAP.html" ,
      " and fill in with the following details to get the number of ions required for the system")
print("Protein mass (kDa):", protein_weight)
print("Solution salt concentration (mM/l):", salt_conc)
print("Net charge of solutes (proton charge units):", int(charge))
print("Number of water molecules:", no_water)
