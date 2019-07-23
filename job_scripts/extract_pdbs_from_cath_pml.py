#!/usr/bin/env python3
'''Extract the PDB files from the pml file from the CATH databse.
The extracted PDBs will be stored in a pdb file.
The PDB files will be cleaned by reading into Rosetta.
'''

import os
import sys

import pyrosetta
from pyrosetta import rosetta


def extract_pdbs_from_cath_pml(pml_file, output_path):
    '''Extract the PDB files from the pml file from the CATH databse.
    The extracted PDBs will be stored in a pdb file. 
    '''

    with open(pml_file, 'r') as f:
        while True:
            line = f.readline()
            if not line: break
         
            # Read one structure

            if line.strip().startswith('cmd.read_pdbstr'):
                pdb_lines = [line.strip()[19:].strip('\\')]
                pdb_id = ''

                while True:  
                    line = f.readline()
                    if line.strip().startswith('"""'):
                        pdb_id = line.strip()[5:12]
                        break

                    pdb_line = line.strip().strip('\\')
                    if len(pdb_line) > 17:
                        pdb_line = pdb_line[0:16] + ' ' + pdb_line[17:] # Remove all altLoc flags
                    
                    pdb_lines.append(pdb_line) # Remove all altLoc flags
           

                try:
                    pose = rosetta.core.pose.Pose()
                    rosetta.core.import_pose.pose_from_pdbstring(pose, '\n'.join(pdb_lines))
                    
                    pose.dump_pdb(os.path.join(output_path, '{0}.pdb.gz'.format(pdb_id)))

                except:
                    continue

def extract_all_pdbs(input_path, output_path, num_jobs, job_id):
    pml_files = [f for f in os.listdir(input_path) if f.endswith('.pml')]
    
    for i, f in enumerate(pml_files):
        if job_id == i % num_jobs:
            os.makedirs(os.path.join(output_path, f[:-4]), exist_ok=True) 
           
            extract_pdbs_from_cath_pml(os.path.join(input_path, f), os.path.join(output_path, f[:-4]))        



if __name__ == '__main__':
    
    data_path = sys.argv[1]
    input_path = sys.argv[2]
    
    num_jobs = 1
    job_id = 0
    if len(sys.argv) > 4:
      num_jobs = int(sys.argv[3])
      job_id = int(sys.argv[4]) - 1


    pyrosetta.init(options='-ignore_unrecognized_res true')

    extract_all_pdbs(input_path, data_path, num_jobs, job_id)
