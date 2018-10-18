#!/usr/bin/env python3
'''Extract the PDB files from the pml file from the CATH databse.
The extracted PDBs will be stored in a pdb file.
The PDB files will be cleaned by reading into Rosetta.

Usage:
    ./extract_pdbs_from_cath_pml.py pml_file
'''

import os
import sys

import pyrosetta
from pyrosetta import rosetta


def extract_pdbs_from_cath_pml(pml_file):
    '''Extract the PDB files from the pml file from the CATH databse.
    The extracted PDBs will be stored in a pdb file. 
    '''
    os.makedirs('pdbs', exist_ok=True) 

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
                    
                    pose.dump_pdb(os.path.join('pdbs', '{0}.pdb.gz'.format(pdb_id)))

                except:
                    continue

if __name__ == '__main__':
    pml_file = sys.argv[1]

    pyrosetta.init(options='-ignore_unrecognized_res true')

    extract_pdbs_from_cath_pml(pml_file) 
