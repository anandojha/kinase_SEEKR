from scipy.interpolate import make_interp_spline
from biopandas.pdb import PandasPdb
import matplotlib.pyplot as plt
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import pickle as pk
import pandas as pd
import numpy as np
import itertools
import warnings
import parmed
import re
import os
warnings.filterwarnings('ignore')
####################################################################################################
# System Parameters
ligand = 'F9J'
pos_ion_1 = 'Na+'
neg_ion_1 = 'Cl-'
distance = "3.0"
file_to_read = "system_nvt_output_last_frame.pdb"
####################################################################################################
pdb_file = "system_final_analyse.pdb"
os.system("rm -rf saved_pdbs")
os.system("mkdir saved_pdbs")
fin = open(file_to_read, "rt")
fout = open("system_V.pdb", "wt")
for line in fin:
    fout.write(line.replace('HETATM', 'ATOM  '))
fin.close()
fout.close()
ppdb = PandasPdb()
ppdb.read_pdb("system_V.pdb")
ppdb.df['ATOM'].loc[(ppdb.df['ATOM'].chain_id == 'A'),'chain_id']= ''
ppdb.to_pdb(path = "system_VI.pdb", records=None, gz=False, append_newline=True)
file1 = open('system_VI.pdb', 'r') 
file2 = open('system_VII.pdb', 'w')  
for line in file1.readlines(): 
    if not (line.startswith('CONECT')): 
        file2.write(line) 
file2.close() 
file1.close()
file1 = open('system_VII.pdb', 'r') 
file2 = open('system_VIII.pdb', 'w')  
for line in file1.readlines(): 
    if not (line.startswith('ENDMDL')): 
        file2.write(line) 
file2.close() 
file1.close()
os.system("pdb4amber -i system_VIII.pdb -o system_IX.pdb --noter")
fin = open("system_IX.pdb", "rt")
fout = open("system_X.pdb", "wt")
for line in fin:
    fout.write(line.replace('HETATM', 'ATOM  '))
fin.close()
fout.close()
os.system("rm -rf *renum* *nonprot* *sslink")

command = "rm -rf system_V.pdb system_VI.pdb system_VII.pdb system_VIII.pdb system_IX.pdb"
os.system(command)

ppdb = PandasPdb()
ppdb.read_pdb("system_X.pdb")
list_new = []
for i in range(len(ppdb.df['ATOM'])):
    list_new.append(i+1)
df_new = pd.DataFrame(list_new) 
df_new.columns = ["atom_number"]
ppdb.df['ATOM']["atom_number"] = df_new["atom_number"]
ppdb.to_pdb(path = pdb_file, records=None, gz=False, append_newline=True)
ppdb0 = PandasPdb()
ppdb0.read_pdb(pdb_file)
df_atoms = ppdb0.df['ATOM']
df_ligand = df_atoms.loc[df_atoms['residue_name'] == ligand]
df_ligand_list = df_ligand["atom_number"]
lig_list = df_ligand_list.values.tolist()
print("Ligand Atom Numbers : ", lig_list)

ppdb1 = PandasPdb().read_pdb(pdb_file)
df = ppdb1.df['ATOM'][['atom_number','residue_number','residue_name']]
lig_index_list = []
for i in lig_list:
    indices = np.where(df["atom_number"] == i)
    indices = list(indices)[0]
    indices = list(indices)
    lig_index_list.append(indices)
lig_index_list = list(itertools.chain.from_iterable(lig_index_list))
print("Ligand Atom Indices : ", lig_index_list)
df_ligand = df_atoms.loc[df_atoms['residue_name'] == ligand]
lig_df= df_ligand[['x_coord', 'y_coord', 'z_coord']]
print(lig_df.head())
ligand_coordinate_list = lig_df.values.tolist()

# Create a dataframe without water and ligand
ppdb2 = PandasPdb().read_pdb(pdb_file)
ppdb2.df['ATOM'] = ppdb2.df['ATOM'][ppdb2.df['ATOM']['residue_name'] != 'WAT']
ppdb2.df['ATOM'] = ppdb2.df['ATOM'][ppdb2.df['ATOM']['residue_name'] != 'HOH']
ppdb2.df['ATOM'] = ppdb2.df['ATOM'][ppdb2.df['ATOM']['residue_name'] != pos_ion_1]
ppdb2.df['ATOM'] = ppdb2.df['ATOM'][ppdb2.df['ATOM']['residue_name'] != neg_ion_1]
ppdb2.df['ATOM'] = ppdb2.df['ATOM'][ppdb2.df['ATOM']['residue_name'] != ligand]
ppdb2.to_pdb(path='saved_pdbs/receptor_dry.pdb', records=None, gz=False, append_newline=True)
ppdb3 = PandasPdb().read_pdb('saved_pdbs/receptor_dry.pdb')
df_atoms_new = ppdb3.df['ATOM']
print(ppdb3.df['ATOM'].shape)
ppdb3.df['ATOM'].head()
rec_list = []
for i in range(len(ligand_coordinate_list)):
    reference_point = ligand_coordinate_list[i]
    ppdb3 = PandasPdb().read_pdb('saved_pdbs/receptor_dry.pdb')
    distances = ppdb3.distance(xyz=reference_point, records=('ATOM',))
    all_within_distance = ppdb3.df['ATOM'][distances < float(distance)]
    receptor_df = all_within_distance["atom_number"]
    receptor_list = receptor_df.values.tolist()
    rec_list.append(receptor_list)
rec_list = list(itertools.chain(*rec_list))
rec_list = set(rec_list)
rec_list = list(rec_list)
rec_list.sort() 
print(rec_list)

ppdb4 = PandasPdb().read_pdb('saved_pdbs/receptor_dry.pdb')
df = ppdb4.df['ATOM'][['atom_number','residue_number','residue_name']]
index_list = []
for i in rec_list:
    indices = np.where(df["atom_number"] == i)
    indices = list(indices)[0]
    indices = list(indices)
    index_list.append(indices)
index_list = list(itertools.chain.from_iterable(index_list))
print(index_list)
df1 = df.iloc[index_list,] 
print(df1.shape)
resid_num = list(df1.residue_number.unique())
print(resid_num)
ppdb5 = PandasPdb().read_pdb('saved_pdbs/receptor_dry.pdb')
df = ppdb5.df['ATOM'][['atom_number','residue_number','residue_name']]
rec_index_list_list = []
for i in resid_num:
    indices = np.where(df["residue_number"] == i)
    indices = list(indices)[0]
    indices = list(indices)
    rec_index_list_list.append(indices)  
rec_index_list = list(itertools.chain.from_iterable(rec_index_list_list))
print("rec_index_list_list")
print(rec_index_list_list)
print("rec_index_list")
print(rec_index_list)
df1 = df.iloc[rec_index_list,] 
print(df1.shape)
df2 = df1['atom_number']
recep_list = df2.values.tolist()
print("recep_list")
print(recep_list)
recep_list_list = []
for i in range(len(rec_index_list_list)):
    add_list = [x + 1 for x in rec_index_list_list[i]]
    recep_list_list.append(add_list)  
print("recep_list_list")
print(recep_list_list)
selected_atoms = []
selected_atoms.extend(lig_list)
selected_atoms.extend(recep_list)
print(selected_atoms)
ppdb6 = PandasPdb().read_pdb(pdb_file)
atom_number_index = []
for i in range(len(ppdb6.df['ATOM'])):
    atom_number_index.append(i+1)
non_selected_region = list(set(atom_number_index).difference(selected_atoms))
print(len(non_selected_region))
print(len(selected_atoms))
print(len(atom_number_index))
print(len(non_selected_region) + len(selected_atoms) == len(atom_number_index))
ppdb7 = PandasPdb().read_pdb(pdb_file)
for i in non_selected_region:
    ppdb7.df['ATOM'] = ppdb7.df['ATOM'][ppdb7.df['ATOM']['atom_number'] != i]
ppdb7.to_pdb(path='saved_pdbs/selected_region.pdb', records=None, gz=False, append_newline=True)
ppdb8 = PandasPdb().read_pdb('saved_pdbs/selected_region.pdb')
df_atoms_new = ppdb8.df['ATOM']
print(ppdb8.df['ATOM'].shape)
# Ligand Pdb 
ligand_atoms = []
ligand_atoms.extend(lig_list)
print(ligand_atoms)
ppdb9 = PandasPdb().read_pdb(pdb_file)
atom_number_index = []
for i in range(len(ppdb9.df['ATOM'])):
    atom_number_index.append(i+1)
non_ligand_selected_region = list(set(atom_number_index).difference(ligand_atoms))
print(len(non_ligand_selected_region))
print(len(ligand_atoms))
print(len(atom_number_index))
print(len(non_ligand_selected_region) + len(ligand_atoms) == len(atom_number_index))
ppdb10 = PandasPdb().read_pdb(pdb_file)
for i in non_ligand_selected_region:
    ppdb10.df['ATOM'] = ppdb10.df['ATOM'][ppdb10.df['ATOM']['atom_number'] != i]
ppdb10.to_pdb(path='saved_pdbs/ligand.pdb', records=None, gz=False, append_newline=True)
ppdb11 = PandasPdb().read_pdb('saved_pdbs/ligand.pdb')
df_lig = ppdb11.df['ATOM']
print(ppdb11.df['ATOM'].shape)
# Ligand Surrounding Pdb 
ligand_sur_atoms = []
ligand_sur_atoms.extend(recep_list)
print(ligand_sur_atoms)
ppdb12 = PandasPdb().read_pdb(pdb_file)
atom_number_index = []
for i in range(len(ppdb12.df['ATOM'])):
    atom_number_index.append(i+1)
non_ligand_sur_selected_region = list(set(atom_number_index).difference(ligand_sur_atoms))
print(len(non_ligand_sur_selected_region))
print(len(ligand_sur_atoms))
print(len(atom_number_index))
print(len(non_ligand_sur_selected_region) + len(ligand_sur_atoms) == len(atom_number_index))
ppdb13 = PandasPdb().read_pdb(pdb_file)
for i in non_ligand_sur_selected_region:
    ppdb13.df['ATOM'] = ppdb13.df['ATOM'][ppdb13.df['ATOM']['atom_number'] != i]
ppdb13.to_pdb(path='saved_pdbs/ligand_sur.pdb', records=None, gz=False, append_newline=True)
ppdb14 = PandasPdb().read_pdb('saved_pdbs/ligand_sur.pdb')
df_lig_sur = ppdb14.df['ATOM']
print(ppdb14.df['ATOM'].shape)
# Create PDB of receptors with alpha carbon atoms surrounding the ligand
ppdb17 = PandasPdb().read_pdb('saved_pdbs/ligand_sur.pdb')
ppdb17.df['ATOM'] = ppdb17.df['ATOM'] [ppdb17.df['ATOM'].atom_name.str.startswith(('CA'))]
print(ppdb17.df['ATOM'].shape)
ppdb17.to_pdb(path='saved_pdbs/ligand_sur_CA.pdb', records=None, gz=False, append_newline=True)
ppdb18 = PandasPdb().read_pdb('saved_pdbs/ligand_sur_CA.pdb')
print(ppdb18.df['ATOM'].shape)
recep_list_CA = ppdb18.df['ATOM']['atom_number'].tolist()
print("Receptor Alpha-Carbon Atom Numbers : ", recep_list_CA)
ppdb19 = PandasPdb().read_pdb(pdb_file)
df = ppdb19.df['ATOM'][['atom_number','residue_number','residue_name']]
recep_list_CA_index = []
for i in recep_list_CA:
    indices = np.where(df["atom_number"] == i)
    indices = list(indices)[0]
    indices = list(indices)
    recep_list_CA_index.append(indices)
recep_list_CA_index = list(itertools.chain.from_iterable(recep_list_CA_index))
print("Receptor Alpha-Carbon Atom Indices : ", recep_list_CA_index)
with open("ligand_receptor_indices.txt", "w") as f:
    f.write("lig_indices : " + "\n")
    f.write(str(lig_index_list) + "\n")
    f.write("rec_indices : " + "\n")
    f.write(str(recep_list_CA_index) + "\n")
####################################################################################################
def get_lig_rec_distance(parmed_struct, positions, lig_atom_list, rec_atom_list):
    parmed_struct.coordinates = positions
    center_of_mass_1 = np.array([[0., 0., 0.]])
    total_mass = 0.0
    for atom_index in lig_atom_list:
        atom_pos = parmed_struct.coordinates[atom_index,:]
        atom_mass = parmed_struct.atoms[atom_index].mass
        center_of_mass_1 += atom_mass * atom_pos
        total_mass += atom_mass
    center_of_mass_1 = center_of_mass_1 / total_mass
    
    center_of_mass_2 = np.array([[0., 0., 0.]])
    total_mass = 0.0
    for atom_index in rec_atom_list:
        atom_pos = parmed_struct.coordinates[atom_index,:]
        atom_mass = parmed_struct.atoms[atom_index].mass
        center_of_mass_2 += atom_mass * atom_pos
        total_mass += atom_mass
    center_of_mass_2 = center_of_mass_2 / total_mass
    distance = np.linalg.norm(center_of_mass_2 - center_of_mass_1)
    return distance
####################################################################################################
lines_indices = open("ligand_receptor_indices.txt",'r').read().splitlines()
ligand_indices = lines_indices[1]
ligand_indices = ligand_indices.strip('][').split(', ')
ligand_indices_list = []
for i in ligand_indices:
    i  = int(i)
    ligand_indices_list.append(i)
receptor_indices = lines_indices[3]
receptor_indices = receptor_indices.strip('][').split(', ')
receptor_indices_list = []
for i in receptor_indices:
    i  = int(i)
    receptor_indices_list.append(i)
print("Length of the receptor alpha carbon list is : ", len(receptor_indices_list))
####################################################################################################
lig_indices = ligand_indices_list
rec_indices = receptor_indices_list
prmtop_filename = "system_TP4EW_I.parm7"
pdb_filename = "system_X.pdb"
####################################################################################################
mypdb = PDBFile(pdb_filename)
parmed_struct = parmed.load_file(pdb_filename, structure = True)
####################################################################################################
with open('rec_lig_dis.out', "w") as f:
    for L in range(15, len(rec_indices) + 1):
        for comb in itertools.combinations(rec_indices, L):
            rec_index_list = list(comb)
            lig_rec_distance = get_lig_rec_distance(parmed_struct = parmed_struct, positions = mypdb.positions, lig_atom_list = lig_indices, rec_atom_list = rec_index_list)
            print(lig_rec_distance, rec_index_list)
            f.write(str(lig_rec_distance) + " " + str(rec_index_list) + "\n")
####################################################################################################
file = open('rec_lig_dis.out', 'r')
lines = file.readlines()
lines = lines[:]
with open('analyse_lig_rec.out', "w") as f:
    for x in lines:
        lig_rec_dist = float(re.findall('\d*\.?\d+', x)[0])
        if lig_rec_dist < 2.51 and lig_rec_dist > 1.99 :
            print(lig_rec_dist)
            f.write(str(lig_rec_dist)  + "\n") 
            rec_list = re.findall('\d*\.?\d+', x)[1:]
            rec_list = [int(x) for x in rec_list]
            f.write(str(rec_list)  + "\n")
            print(rec_list)   
#################################################################################################### 

