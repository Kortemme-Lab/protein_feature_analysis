#!/usr/bin/env python3
'''Relax the structures using pyrosetta
'''

import sys
import os
from datetime import timedelta

import numpy as np
import pandas as pd
from flufl.lock import Lock

import pyrosetta
from pyrosetta import rosetta

def append_to_csv(data_list, file_name):
    '''Append a list of tuples to a csv file. This function is thread safe.'''
    
    # Make a lock

    lock_name = os.path.join(file_name + '.lock')
    lock = Lock(lock_name)
    lock.lifetime = timedelta(minutes=10)

    # Write to the result file with a lock
    
    with lock:
      with open(file_name, 'a+') as f:
        for t in data_list:
            f.write(','.join([str(x) for x in t]) + '\n')


def get_superimpose_transformation(P1, P2):
    '''Get the superimpose transformation that transfoms a list of
    points P1 to another list of points P2.'''
    if len(P1) != len(P2):
        raise Exception("Sets to be superimposed must have same number of points.")

    com1 = np.mean(P1, axis=0)
    com2 = np.mean(P2, axis=0)

    R = np.dot(np.transpose(np.array(P1) - com1), np.array(P2) - com2)
    V, S, W = np.linalg.svd(R)

    if (np.linalg.det(V) * np.linalg.det(W)) < 0.0:
        V[:, -1] = -V[:, -1]

    M = np.transpose(np.array(np.dot(V, W)))

    return M, com2 - np.dot(M, com1)

def align_points_and_calculate_RMSD(xyzs1, xyzs2):
    '''Align the points and calculate the RMSD after alignment.'''
    def xyzV_to_np_array(xyz):
        return np.array([xyz.x, xyz.y, xyz.z])

    def RMSD(P1, P2):
        diff = [P1[i] - P2[i] for i in range(len(P1))]
        return np.sqrt(sum(np.dot(d, d) for d in diff) / len(P1))

    P1 = [xyzV_to_np_array(xyz) for xyz in xyzs1]
    P2 = [xyzV_to_np_array(xyz) for xyz in xyzs2]

    M, v = get_superimpose_transformation(P1, P2)

    P1_aligned = [np.dot(M, p) + v for p in P1]

    return RMSD(P1_aligned, P2)

def relax_one_structure(pose, relax_bb=True):
    '''Relax one structure'''
    # MinPack
    
    task_factory = rosetta.core.pack.task.TaskFactory()
    task_factory.push_back(rosetta.core.pack.task.operation.RestrictToRepacking())
    
    min_pack = rosetta.protocols.minimization_packing.MinPackMover()
    min_pack.task_factory(task_factory)
    min_pack.apply(pose)

    # Minimize

    mm = rosetta.core.kinematics.MoveMap()
    mm.set_bb(relax_bb)
    mm.set_chi(True)
    mm.set_jump(relax_bb)

    min_mover = rosetta.protocols.minimization_packing.MinMover()
    min_mover.set_movemap(mm)
    min_mover.apply(pose)

def relax_structures(data_path, input_path, num_jobs, job_id):
    '''Relax the structures in a input path and save the 
    relaxed structures to the data_path. 
    Also write the backbone RMSD difference between the original and
    relaxed structures.
    '''
    pyrosetta.init('-ignore_unrecognized_res true -packing:use_input_sc true')

    for i, f in enumerate(os.listdir(input_path)):
        if i % num_jobs != job_id:
            continue

        if f.endswith('.pdb') or f.endswith('.pdb.gz'):
            try:
            
                pose = rosetta.core.pose.Pose()
                rosetta.core.import_pose.pose_from_file(pose, os.path.join(input_path, f))

                bb_atoms = ['N', 'CA', 'C', 'O']
                xyzs_before = [pose.residue(i).xyz(a) for i in range(1, pose.size() + 1) for a in range(1, pose.residue(i).n_mainchain_atoms())] 
                
                relax_one_structure(pose)
                xyzs_after = [pose.residue(i).xyz(a) for i in range(1, pose.size() + 1) for a in range(1, pose.residue(i).n_mainchain_atoms())] 
               
                bb_rmsds = align_points_and_calculate_RMSD(xyzs_before, xyzs_after)

                data_list = [(f.split('.')[0], bb_rmsds)]
                append_to_csv(data_list, os.path.join(data_path, 'summary.csv'))

                pose.dump_pdb(os.path.join(data_path, f.split('.')[0] + '.pdb.gz'))

            except:
                continue

if __name__ == '__main__':
    
    data_path = sys.argv[1]
    input_path = sys.argv[2]
    
    num_jobs = 1
    job_id = 0
    if len(sys.argv) > 4:
      num_jobs = int(sys.argv[3])
      job_id = int(sys.argv[4]) - 1
    
    relax_structures(data_path, input_path, num_jobs, job_id)
