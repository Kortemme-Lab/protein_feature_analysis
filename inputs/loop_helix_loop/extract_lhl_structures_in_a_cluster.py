#!/usr/bin/env python3
'''Extract LHL structures in a cluster
Save the structures into a directory called clustered_lhl_structure.
Also save the information of insertion points.

Usage:
    ./extract_lhl_structures_in_a_cluster.py pdbs_path cluster_file
'''

import os
import sys
import json

import pyrosetta
from pyrosetta import rosetta


def slice_peptide(pose, start, stop):
    '''Slice a peptide from a pose into a new pose.'''
    seqposes = rosetta.utility.vector1_unsigned_long()
    for seqpos in range(start, stop + 1):
        seqposes.append(seqpos)
    
    pp_pose = rosetta.core.pose.Pose()
    rosetta.core.pose.pdbslice(pp_pose, pose, seqposes)
   
    return pp_pose

def find_ss_region_of_position_pose(pose, seqpos):
    '''Find the start and stop of the secondary structure
    element that a seqpos is within.
    '''
    dssp_str = rosetta.core.scoring.dssp.Dssp(pose).get_dssp_secstruct()
   
    start = seqpos
    while start > 1 and dssp_str[start - 1] == dssp_str[start - 2]:
        start -= 1
    
    stop = seqpos 
    while stop < pose.size() and dssp_str[stop - 1] == dssp_str[stop]:
        stop += 1

    return start, stop

def extract_lhl_structures_in_a_cluster(pdbs_path, cluster_file):
    '''Get LHL distributions'''
    # Load the cluster

    with open(cluster_file, 'r') as f:
        cluster = json.load(f)

    # Dump the lhl units in the cluster

    #print(len(cluster))

    flanking_residues = []
    os.makedirs('clustered_lhl_structure', exist_ok=True)

    for i, lhl_unit in enumerate(cluster):
        pose = rosetta.core.import_pose.pose_from_file(os.path.join(pdbs_path, lhl_unit['pdb_file'])) 

        pre_ss_start = find_ss_region_of_position_pose(pose, lhl_unit['start'] - 1)[0]       
        post_ss_stop = find_ss_region_of_position_pose(pose, lhl_unit['stop'] + 1)[1]    

        pp_pose = slice_peptide(pose, pre_ss_start, post_ss_stop)
        pp_pose.dump_pdb(os.path.join('clustered_lhl_structure', 'model_{0}.pdb.gz'.format(i)))

        with open(os.path.join('clustered_lhl_structure', 'insertion_points_{0}.json'.format(i)), 'w') as f:
            json.dump([{'start':lhl_unit['start'] - pre_ss_start + 1, 'stop':lhl_unit['stop'] - pre_ss_start + 1}], f)

        for j in range(pre_ss_start, lhl_unit['start']):
            flanking_residues.append((i, pose.pdb_info().pose2pdb(j).split(' ')[0]))

        for j in range(lhl_unit['stop'] + 1, post_ss_stop + 1):
            flanking_residues.append((i, pose.pdb_info().pose2pdb(j).split(' ')[0]))

    # Print the PyMol selection command for the flanking residues

    pymol_str = 'sele flanking_residues,'

    for i, j in flanking_residues:
        pymol_str += ' (model_{0} and res {1})'.format(i, j)

    print(pymol_str)

if __name__ == '__main__':
    pdbs_path = sys.argv[1]
    cluster_file = sys.argv[2]

    pyrosetta.init(options='-ignore_unrecognized_res true')

    extract_lhl_structures_in_a_cluster(pdbs_path, cluster_file)
