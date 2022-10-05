from biopandas.pdb import PandasPdb
from simtk.openmm import app
from sys import stdout
from time import time
from prepare import *
import pandas as pd
import numpy as np
import parmed
import simtk
import time
import os
import re

lines_indices = open("ligand_receptor_indices.txt",'r').read().splitlines()
ligand_indices = lines_indices[1]
ligand_indices = ligand_indices.strip('][').split(', ')
ligand_indices_list = []
for i in ligand_indices:
    i  = int(i)
    ligand_indices_list.append(i)
file = open('analyse_lig_rec.out', 'r')
lines = file.readlines()
dist_lines  = lines[ : : 2]
rec_list_lines  = lines[1 : : 2]
dist_list = []
for x in dist_lines:
    lig_rec_dist = float(re.findall('\d*\.?\d+', x)[0])
    dist_list.append(lig_rec_dist)
rec_atom_list = []
for x in rec_list_lines:
    rec_list = re.findall('\d*\.?\d+', x)[:]
    rec_list = [int(x) for x in rec_list]
    rec_atom_list.append(rec_list)
for i in range(len(dist_list)):
    if dist_list[i] == min (dist_list):
        rec_index_selected_list = rec_atom_list[i]
lig_indices = ligand_indices_list
rec_indices = rec_index_selected_list

restrainToBoundState(
    rec_indices = rec_indices, 
    lig_indices = lig_indices,
    num_steps=10000000,
    steps_per_energy_update=5000,
    steps_per_trajectory_update=50000,
    time_step=0.002,
    spring_constant=3000.0,
    target_radius=1.00,
    temperature=300.00,
    cuda_index=0,
    nonbonded_cutoff=0.9,
    raw_prmtop_filename="system_TP4EW_I.parm7",
    raw_inpcrd_filename="system_TP4EW_I.inpcrd",
    prmtop_filename="hostguest.parm7",
    inpcrd_filename="hostguest_raw2.rst7",
    inpcrd_final="hostguest.rst7",
    input_file="system_X.pdb",
    trajectory_filename="hostguest_restrain_bound_state.pdb",
    output_file="hostguest_at2.50.pdb",
)

pullOutLigandSMD(
    rec_indices = rec_indices, 
    lig_indices = lig_indices,
    target_radii= [2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00],
    prmtop_filename="hostguest.parm7",
    inpcrd_filename="hostguest.rst7",
    pdb_filename="hostguest_at2.50.pdb",
    temperature=300.00,
    spring_constant=50000.0,
    nonbonded_cutoff=0.9,
    trajectory_filename="hostguest_out.pdb",
    trajectory_frequency=100000,
    total_num_steps=500000000,
    num_windows=100,
    show_state_output=False,
    state_filename="state.xml",
    basename="hostguest_at",
    cuda_device_index=0,
)
