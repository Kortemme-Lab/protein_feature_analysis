#!/usr/bin/env python3
'''Find the LHL units of pdb files in a given path.
Dump the LHL definition files into a lhl_info path.

Usage:
    ./find_lhl_units_for_pdbs.py pdbs_path
'''

import os
import sys
import json

import pyrosetta
from pyrosetta import rosetta


def find_lhl_units_for_one_pdb(pdb_file):
    '''Find LHL units for one pdb file.
    Return a list of the LHL units information.
    '''
    pose = rosetta.core.import_pose.pose_from_file(pdb_file)
    dssp_str = rosetta.core.scoring.dssp.Dssp(pose).get_dssp_secstruct()

    lhl_info = []

    # Aggregate secondary structures

    ss_structs = []

    start = 0

    while start < len(dssp_str):
        ss = dssp_str[start]
        stop = start
       
        while stop < len(dssp_str) and ss == dssp_str[stop]:
            stop += 1

        ss_structs.append((start + 1, stop, ss))
   
        start = stop

    # Find all LHL units

    for i in range(1, len(ss_structs) - 3):
        if ss_structs[i][2] + ss_structs[i + 1][2] + ss_structs[i + 2][2] == 'LHL':
            start = ss_structs[i][0]
            stop = ss_structs[i + 2][1]
           
            pre_anchor_ca = pose.residue(start - 1).xyz('CA')
            post_anchor_ca = pose.residue(stop + 1).xyz('CA')

            lhl_info.append({'start' : start,
                             'stop' : stop,
                             'H_start' : ss_structs[i + 1][0],
                             'H_stop' : ss_structs[i + 1][1],
                             'ss_pre' : ss_structs[i - 1][2],
                             'ss_post' : ss_structs[i + 3][2],
                             'pdb_file' : os.path.basename(pdb_file),
                            
                             'pre_anchor_ca' : [pre_anchor_ca.x, pre_anchor_ca.y, pre_anchor_ca.z], 
                             'post_anchor_ca' : [post_anchor_ca.x, post_anchor_ca.y, post_anchor_ca.z], 
                            })

    return lhl_info

def find_lhl_units_for_pdbs(pdbs_path):
    '''Find the LHL units of pdb files in a given path.
    Dump the LHL definition files into a lhl_info path.
    '''
    os.makedirs('lhl_info', exist_ok=True)

    for pdb_file in os.listdir(pdbs_path):
        lhl_info =  find_lhl_units_for_one_pdb(os.path.join(pdbs_path, pdb_file))

        with open(os.path.join('lhl_info', '{0}.json'.format(pdb_file.split('.')[0])), 'w') as f:
            json.dump(lhl_info, f)

if __name__ == '__main__':
    pdbs_path = sys.argv[1]

    pyrosetta.init(options='-ignore_unrecognized_res true')

    find_lhl_units_for_pdbs(pdbs_path)
