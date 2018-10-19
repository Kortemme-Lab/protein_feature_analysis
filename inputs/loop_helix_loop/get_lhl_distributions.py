#!/usr/bin/env python3
'''Get distributions for LHL units.
Dump the distributions of various features to text files.

Usage:
    ./get_lhl_distributions.py pdbs_path lhl_info_path edges_file
'''

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

def get_helix_direction(pose, helix_start, helix_stop):
    '''Get the helix direction.
    The direction is defined as the average of the C-O vectors.
    '''
    c_o_vectors = [pose.residue(i).xyz('O') - pose.residue(i).xyz('C')
                    for i in range(helix_start, helix_stop + 1)] 

    sum_vecs = c_o_vectors[0]

    for i in range(1, len(c_o_vectors)):
        sum_vecs += c_o_vectors[i]

    return sum_vecs.normalized()

def get_lhl_lengths(lhl_infos):
    '''Get the distribution of LHL lengths.'''
    return [lhl['stop'] - lhl['start'] + 1 for lhl in lhl_infos]

def get_front_loop_lengths(lhl_infos):
    '''Get the distribution of front loop lengths of LHL units.'''
    return [lhl['H_start'] - lhl['start'] for lhl in lhl_infos]

def get_back_loop_lengths(lhl_infos):
    '''Get the distribution of back loop lengths of LHL units.'''
    return [lhl['stop'] - lhl['H_stop'] for lhl in lhl_infos]

def get_lhl_pair_length_diffs(lhl_infos, edges):
    '''Get the length difference between pairs of LHL units.'''
    length_diffs = []

    for i, j in edges:
        length1 = lhl_infos[i]['stop'] - lhl_infos[i]['start'] + 1
        length2 = lhl_infos[j]['stop'] - lhl_infos[j]['start'] + 1

        length_diffs.append(np.absolute(length1 - length2))

    return length_diffs

def get_lhl_pair_rmsds(poses_map, lhl_infos, edges):
    '''Get the backbone RMSDs between pairs of LHL units'''
    rmsds = []

    for i, j in edges:
        length1 = lhl_infos[i]['stop'] - lhl_infos[i]['start'] + 1
        length2 = lhl_infos[j]['stop'] - lhl_infos[j]['start'] + 1

        len_comp = min(length1, length2)

        residues1 = [lhl_infos[i]['start'] + k for k in range(len_comp)]
        residues2 = [lhl_infos[j]['start'] + k for k in range(len_comp)]

        rmsds.append(calc_backbone_RMSD(poses_map[lhl_infos[i]['pdb_file']], residues1,
            poses_map[lhl_infos[j]['pdb_file']], residues2))

    return rmsds

def get_lhl_pair_helicies_angles(poses_map, lhl_infos, edges):
    '''Get the angles between helices of pairs of LHL units'''
    angles = []

    for i, j in edges:
        helix_direction1 = get_helix_direction(poses_map[lhl_infos[i]['pdb_file']], lhl_infos[i]['H_start'], lhl_infos[i]['H_stop'])
        helix_direction2 = get_helix_direction(poses_map[lhl_infos[j]['pdb_file']], lhl_infos[j]['H_start'], lhl_infos[j]['H_stop'])
      
        cos_angle = helix_direction1.dot(helix_direction2)

        angles.append(180 / np.pi * np.arccos(cos_angle))

    return angles

def dump_distribution(data, data_name):
    '''Dump a distribution to a text file'''

    with open('{0}.txt'.format(data_name), 'w') as f:
        for d in data:
            f.write('{0}\n'.format(d))

def get_lhl_distributions(pdbs_path, lhl_info_path, edges_file):
    '''Get LHL distributions'''
    # Load the pdbs

    poses_map = {}

    for pdb_file in os.listdir(pdbs_path):
        poses_map[pdb_file] = rosetta.core.import_pose.pose_from_file(os.path.join(pdbs_path, pdb_file))

    # Load the lhl_infos
    
    lhl_infos = []

    for lhl_info_file in os.listdir(lhl_info_path):

        with open(os.path.join(lhl_info_path, lhl_info_file), 'r') as f:
            lhl_info = json.load(f)

            lhl_infos += lhl_info

    # Load the edges

    with open(edges_file, 'r') as f:
        edges = json.load(f)

    # Calcualte and dump the distributions

    lhl_lengths = get_lhl_lengths(lhl_infos)
    dump_distribution(lhl_lengths, 'lhl_lengths')

    front_loop_lengths = get_front_loop_lengths(lhl_infos)
    dump_distribution(front_loop_lengths, 'front_loop_lengths')

    back_loop_lengths = get_back_loop_lengths(lhl_infos)
    dump_distribution(back_loop_lengths, 'back_loop_lengths')

    lhl_pair_length_diffs = get_lhl_pair_length_diffs(lhl_infos, edges)
    dump_distribution(lhl_pair_length_diffs, 'lhl_pair_length_diffs')

    lhl_pair_rmsds = get_lhl_pair_rmsds(poses_map, lhl_infos, edges)
    dump_distribution(lhl_pair_rmsds, 'lhl_pair_rmsds')

    lhl_pair_helices_angles = get_lhl_pair_helicies_angles(poses_map, lhl_infos, edges)
    dump_distribution(lhl_pair_helices_angles, 'lhl_pair_helices_angles')

if __name__ == '__main__':
    pdbs_path = sys.argv[1]
    lhl_info_path = sys.argv[2]
    edges_file = sys.argv[3]

    pyrosetta.init(options='-ignore_unrecognized_res true')

    get_lhl_distributions(pdbs_path, lhl_info_path, edges_file)
