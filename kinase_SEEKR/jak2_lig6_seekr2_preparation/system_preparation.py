import os
# Antechamber Ligand Treatment
os.system("antechamber -i ligand.mol2 -fi mol2 -o ligand.in -fo prepi -c bcc -nc 0")
os.system("parmchk2 -i ligand.in -o ligand.frcmod -f prepi -a Y")

# Deleting unneeded files
os.system("rm -rf ANTECHAMBER*")
os.system("rm -rf leap.log")
os.system("rm -rf sqm*")
os.system("rm -rf ATOMTYPE.INF PREP.INF")
os.system("mv NEWPDB.PDB ligand_antechambered.pdb")

# Save the tleap script to file that checks for the initial charge of the system

with open('initial_charge.leap', 'w') as f:
    f.write('''
loadamberparams Y2P.frcmod
loadOff Y2P.off
loadamberprep ligand.in
loadamberparams ligand.frcmod
source leaprc.protein.ff14SB
source leaprc.gaff
source leaprc.water.tip4pew
loadAmberParams frcmod.ionsjc_tip4pew
loadAmberParams frcmod.tip4pew
pdb = loadpdb system.pdb
charge pdb
quit
''')

# Check initial charge 
os.system("tleap -f initial_charge.leap")

# Save the tleap script to file that solvates the system
with open('input_TP4EW.leap', 'w') as f:
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
solvateOct pdb TIP4PEWBOX 10
charge pdb
addions2 pdb Na+ 27
addions2 pdb Cl- 22
charge pdb
saveamberparm pdb system_TP4EW.prmtop system_TP4EW.inpcrd
saveamberparm pdb system_TP4EW.parm7 system_TP4EW.rst7
savepdb pdb system_TP4EW.pdb
quit
''')

# Solvate the system
os.system("tleap -f input_TP4EW.leap")

# Save the tleap script to file that checks for the final charge of the system

with open('final_charge.leap', 'w') as f:
    f.write('''
loadamberparams Y2P.frcmod
loadOff Y2P.off
loadamberprep ligand.in
loadamberparams ligand.frcmod
source leaprc.protein.ff14SB
source leaprc.gaff
source leaprc.water.tip4pew
loadAmberParams frcmod.ionsjc_tip4pew
loadAmberParams frcmod.tip4pew
pdb = loadpdb system_TP4EW.pdb
charge pdb
quit
''')

# Check initial charge 
os.system("tleap -f final_charge.leap")

os.system("rm -rf leap.log initial_charge.leap final_charge.leap input_TP4EW.leap")
