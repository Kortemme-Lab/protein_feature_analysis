#!/usr/bin/env python3

import os
import sys
import json

import numpy as np

import pyrosetta
from pyrosetta import rosetta


def xyzV_to_np_array(xyz):
    return np.array([xyz.x, xyz.y, xyz.z])

def get_backbone_points(pose, residues):
    '''Get backbone points for residues in a pose.'''
    points = []

    for res in residues:
        for atom in ['N', 'CA', 'C']:
            points.append(xyzV_to_np_array(pose.residue(res).xyz(atom)))

    return points

def calc_backbone_RMSD(pose1, residues1, pose2, residues2):
    '''Calculate backbone RMSD between two poses for specific positions.'''
    assert(len(residues1) == len(residues2))

    def RMSD(points1, poinsts2):
        '''Calcualte RMSD between two lists of numpy points.'''
        diff = [points1[i] - poinsts2[i] for i in range(len(points1))]
        return np.sqrt(sum(np.dot(d, d) for d in diff) / len(diff))

    points1 = get_backbone_points(pose1, residues1)
    points2 = get_backbone_points(pose2, residues2)

    return RMSD(points1, points2)

def rmsd_between_two_helices(pose1, helix_start1, helix_stop1, pose2, helix_start2, helix_stop2):
    '''Calculate RMSD between two helices'''
    len_comp = min(helix_stop1 - helix_start1, helix_stop2 - helix_start2)

    h_mid_start1 = (helix_start1 + helix_stop1 - len_comp) // 2
    h_mid_start2 = (helix_start2 + helix_stop2 - len_comp) // 2

    maped_residues_1 = [h_mid_start1 + l for l in range(len_comp)]
    maped_residues_2 = [h_mid_start2 + l for l in range(len_comp)]

    return calc_backbone_RMSD(pose1, maped_residues_1, pose2, maped_residues_2)

def closest_rmsds_between_lhl_units(pose1, lhl_units1, pose2, lhl_units2):
    '''Calculate RMSDs between lhl units in two structures.
    Return the RMSDs between corresponding lhl units defined
    by the closest RMSD.
    '''
    if len(lhl_units2) < len(lhl_units1):
        query_pose = pose2
        query_lhls = lhl_units2
        target_pose = pose1
        target_lhls = lhl_units1
    else:
        query_pose = pose1
        query_lhls = lhl_units1
        target_pose = pose2
        target_lhls = lhl_units2

    closest_rmsds = []

    for q_lhl in query_lhls:
        rmsds = []
        for t_lhl in target_lhls:
            rmsds.append(rmsd_between_two_helices(query_pose, q_lhl['H_start'], q_lhl['H_stop'], target_pose, t_lhl['H_start'], t_lhl['H_stop']))

        closest_rmsds.append(min(rmsds))

    return closest_rmsds

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
           
            lhl_info.append({'start' : start,
                             'stop' : stop,
                             'H_start' : ss_structs[i + 1][0],
                             'H_stop' : ss_structs[i + 1][1],
                             'ss_pre' : ss_structs[i - 1][2],
                             'ss_post' : ss_structs[i + 3][2],
                            })

    return lhl_info

def analyze_one_superfamily(pdb_path):
    '''Return a dictionary of the statistics.'''
    lhl_stat = {}

    pdb_files = [os.path.join(pdb_path, f) for f in os.listdir(pdb_path) if f.endswith('.pdb.gz')]
    lhl_units = [find_lhl_units_for_one_pdb(f) for f in pdb_files]

    # Get the median number of LHL units

    n_lhls = [len(lhls) for lhls in lhl_units]
    lhl_stat['median_number_of_lhl_units'] = np.median(n_lhls)

    # Get the lengths of LHL units

    lhl_stat['lhl_lengths'] = [lhl['stop'] - lhl['start'] + 1 for lhls in lhl_units for lhl in lhls]
    lhl_stat['loop_lengths'] = []
        
    for lhls in lhl_units:
        for lhl in lhls:
            lhl_stat['loop_lengths'].append(lhl['H_start'] - lhl['start'])
            lhl_stat['loop_lengths'].append(lhl['stop'] - lhl['H_stop'])

    # Get the helix rmsd between LHL units

    helix_rmsds = []

    for i in range(len(pdb_files)):
        pose1 = rosetta.core.import_pose.pose_from_file(pdb_files[i])
        
        for j in range(i + 1, len(pdb_files)):
            pose2 = rosetta.core.import_pose.pose_from_file(pdb_files[j])

            helix_rmsds += closest_rmsds_between_lhl_units(pose1, lhl_units[i], pose2, lhl_units[j])

    lhl_stat['helix_rmsds'] = helix_rmsds

    return lhl_stat

def analyze_all_superfamilies(output_path, cath_pdb_path, num_jobs, job_id):
    sfs = [f for f in os.listdir(cath_pdb_path) if f.startswith('cath')]

    for i, sf in enumerate(sfs):
        if job_id == i % num_jobs:
            lhl_stat = analyze_one_superfamily(os.path.join(cath_pdb_path, sf))

            with open(os.path.join(output_path, sf + '.json'), 'w') as f:
                json.dump(lhl_stat, f)

if __name__ == '__main__':
    
    data_path = sys.argv[1]
    cath_pdb_path = sys.argv[2]
    
    num_jobs = 1
    job_id = 0
    if len(sys.argv) > 4:
      num_jobs = int(sys.argv[3])
      job_id = int(sys.argv[4]) - 1


    pyrosetta.init(options='-ignore_unrecognized_res true')

    analyze_all_superfamilies(data_path, cath_pdb_path, num_jobs, job_id)
