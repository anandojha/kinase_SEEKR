from biopandas.pdb import PandasPdb
from simtk.openmm import app
from sys import stdout
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
lig_indices_from_0 = []
for i in range(len(ligand_indices_list)):
    lig_indices_from_0.append(i)
lig_indices = ligand_indices_list
rec_indices = rec_index_selected_list
ligand_atom_list = lig_indices_from_0

get_host_guest_pqr(
    prmtop_file="system_TP4EW_I.parm7",
    init_pdb="system_X.pdb",
    ligand_pqr="hostguest_ligand.pqr",
    receptor_pqr="hostguest_receptor.pqr",
    lig_resid="F9J",
)

create_prepare_openmmvt_xml(rec_indices = rec_indices,
                            lig_indices = lig_indices,
                            ligand_atom_list = ligand_atom_list,
                            temperature = 300.00,
                            pressure = 1.0,
                            stepsperanchor = 100000000,
                            outputfrequency = 10000000,
                            parmfile = "hostguest.parm7",
                            rstfile = "hostguest.rst7",
                            receptorpqr = "hostguest_receptor.pqr",
                            ligandpqr = "hostguest_ligand.pqr",
                            filename = "prepare.xml",
)

# Execute to make simulation directory
os.system("python /home/aaojha/seekr2/seekr2/prepare.py prepare.xml --skip_checks")
os.system("rm -rf system_seekr2_files")
os.system("mkdir system_seekr2_files")
os.system("mv analyse_lig_rec.out rec_lig_dis.out execute.out prepare_SMD.out ligand_receptor_indices.txt prepare.xml hostguest_raw2.rst7 hostguest_receptor.pqr hostguest_restrain_bound_state.pdb hostguest.rst7 hostguest_ligand.pqr state.xml hostguest_out.pdb system_final_analyse.pdb hostguest.parm7 system_nvt_output_last_frame.pdb system_TP4EW_I.inpcrd system_TP4EW_I.parm7 system_X.pdb system_seekr2_files")
os.system("mv *hostguest_at* system_seekr2_files")
os.system("mv annealing_simulation_box_vectors.pkl annealing.state equilibration.out ligand_antechambered.pdb ligand.frcmod ligand.in ligand.mol2 ligand.pdb npt_simulation_box_vectors.pkl nvt_simulation_box_vectors.pkl system_annealing_output_last_frame.pdb system_annealing_output.pdb  system_I.pdb  system_II_A.pdb system_II.pdb system_III.pdb system_IV.pdb system_npt_output_last_frame.pdb system_npt_output.pdb system_nvt_output.pdb system.pdb system_TP4EW.inpcrd system_TP4EW_I.pdb system_TP4EW.parm7 system_TP4EW.pdb system_TP4EW.prmtop system_TP4EW.rst7 init_system_5.pdb Y2P.frcmod Y2P.off system_seekr2_files")

current_pwd = os.getcwd()
os.chdir(current_pwd + "/system_seekr2")
lines_file = open("model.xml", "r")
lines = lines_file.readlines()
for i in range(len(lines)):
    if "<num_production_steps" in lines[i]:
        to_begin = int(i)
    if "<restart_checkpoint_interval" in lines[i]:
        to_end = int(i)
lines_I = lines[: to_begin]
lines_II = lines[to_end + 1:]
with open("model_test.xml", "w") as file:
    for i in lines_I:
        file.write(i)
    file.write ("        <num_production_steps type=" + '"' + "int" + '"' + ">10000</num_production_steps>"  + "\n")
    file.write ("        <energy_reporter_interval type=" + '"' + "int" + '"' + ">10000</energy_reporter_interval>"  + "\n")
    file.write ("        <trajectory_reporter_interval type=" + '"' + "int" + '"' + ">10000</trajectory_reporter_interval>"  + "\n")
    file.write ("        <restart_checkpoint_interval type=" + '"' + "int" + '"' + ">10000</restart_checkpoint_interval>"  + "\n")
    for i in lines_II:
        file.write(i)

# Run Test SEEKR2 calculations
anchor_list = []
for i in os.listdir():
    if i.startswith("anchor_"):
        anchor_list.append(i)
for i in range(len(anchor_list)):
    print("Running anchor" + " " + str(i))
    command = "python /home/aaojha/seekr2/seekr2/run.py" + " " + str(i) + " " + "model_test.xml" + " -f" + " > " + str(i) + ".out"
    os.system(command)
    file_to_analyze = str(i) + ".out"
    file_size = os.stat(file_to_analyze).st_size
    if file_size == 0:
        print("Anchor " + str(i) +  " is facing issues need to be fixed")
    else:
        print("All good with this anchor")
    command = "rm -rf " + str(i) + ".out"
    os.system(command)
command = "rm -rf model_test.xml"
os.system(command)
os.chdir (current_pwd)
command = "rm -rf __pycache__"
os.system(command)

print("We are ready to run SEEKR2 calculations.")

#Removing the production steps files in the anchors
current_pwd = os.getcwd()
os.chdir(current_pwd + "/system_seekr2")
anchor_list = []
for i in os.listdir():
    if i.startswith("anchor_"):
        anchor_list.append(i)
for i in anchor_list:
    os.chdir(current_pwd + "/" + "system_seekr2" + "/" + i + "/" + "prod")
    os.system("rm -rf *")
    os.chdir(current_pwd) 



