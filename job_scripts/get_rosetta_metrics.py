#!/usr/bin/env python3
'''Get rosetta related metrics
'''

import sys
import os
from datetime import timedelta
import json

import numpy as np
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

def get_score_data(pose, structure_name):
    '''Get the list of tuples for score.'''
    sfxn = rosetta.core.scoring.get_score_function()
    score = sfxn(pose)
    return [(structure_name, score)]

def get_residue_average_score_data(pose, structure_name):
    '''Get the list of tuples for the residue average score.'''
    sfxn = rosetta.core.scoring.get_score_function()
    score = sfxn(pose)
    return [(structure_name, score / pose.size())]

def get_residue_score_data(pose, structure_name):
    '''Get the list of tuple for scores for each residue'''
    sfxn = rosetta.core.scoring.get_score_function()
    sfxn(pose)
   
    data = []

    for i in range(1, pose.size() + 1):
        data.append(('{0}_{1}_{2}'.format(structure_name, i, pose.residue(i).name3()),
            pose.energies().residue_total_energy(i)))

    return data

def get_num_buried_unsatisfied_hbonds(pose):
    '''Get the number of buried unsatisfied hbonds.'''
    bupc = rosetta.protocols.simple_pose_metric_calculators.BuriedUnsatisfiedPolarsCalculator(
            'default', 'default')

    sfxn = rosetta.core.scoring.get_score_function()
    sfxn(pose)
    
    buhs_for_each_res = json.loads(bupc.get('residue_bur_unsat_polars', pose))

    return sum(buhs_for_each_res)

def get_num_buried_unsatisfied_hbonds_data(pose, structure_name):
    '''Get the list of tuples for the number of buried unsatisfied hbonds.'''
    return [(structure_name, get_num_buried_unsatisfied_hbonds(pose))]

def get_average_num_buried_unsatisfied_hbonds_data(pose, structure_name):
    '''Get the list of tuples for the average number of buried unsatisfied hbonds.'''
    return [(structure_name, get_num_buried_unsatisfied_hbonds(pose) / pose.size())]

def get_num_over_saturated_hbond_acceptors(pose):
    '''Get the number of over saturated hbond acceptors for'''
    oshaf = rosetta.protocols.cyclic_peptide.OversaturatedHbondAcceptorFilter()
    oshaf.set_consider_mainchain_only(False)
    return oshaf.report_sm(pose)

def get_num_over_saturated_hbond_acceptors_data(pose, structure_name):
    '''Get the list of tuples for the number of over saturated hbond acceptors.'''
    return [(structure_name, get_num_over_saturated_hbond_acceptors(pose))]

def get_average_num_over_saturated_hbond_acceptors_data(pose, structure_name):
    '''Get the list of tuples for the average number of over saturated hbond acceptors.'''
    return [(structure_name, get_num_over_saturated_hbond_acceptors(pose) / pose.size())]

def get_hydrophobic_sasa(pose):
    '''Calculate the hydrophobic sasa'''
    rsd_sasa = pyrosetta.rosetta.utility.vector1_double()
    rsd_hydrophobic_sasa = pyrosetta.rosetta.utility.vector1_double()
    rosetta.core.scoring.calc_per_res_hydrophobic_sasa(pose, rsd_sasa, rsd_hydrophobic_sasa, 1.4) #The last arguement is the probe radius

    return sum(rsd_hydrophobic_sasa)

def get_hydrophobic_sasa_data(pose, structure_name):
    '''Get the list of tuples for the hydrophobic sasa data.'''
    return [(structure_name, get_hydrophobic_sasa(pose))]

def get_average_hydrophobic_sasa_data(pose, structure_name):
    '''Get the list of tuples for the average hydrophobic sasa data.'''
    return [(structure_name, get_hydrophobic_sasa(pose) / pose.size())]

def get_holes_score(pose):
    '''Get the holes score for a list of residues.'''
    dalphaball = 'dependencies/dependencies/DAlpahBall/DAlphaBall.gcc'
    rosetta.basic.options.set_file_option('holes:dalphaball', dalphaball)
   
    # Make a lock

    lock_name = os.path.join('hole_score.lock')
    lock = Lock(lock_name)
    lock.lifetime = timedelta(minutes=3)

    # Write to the result file with a lock
    
    with lock:
        hf = rosetta.protocols.simple_filters.HolesFilter()
        return hf.compute(pose)

def get_holes_score_data(pose, structure_name):
    '''Get the list of tuples for the holes score'''
    return [(structure_name, get_holes_score(pose))]

def get_metrics_for_one_pose(pose, structure_name):
    '''Get rosetta metrics for one pose.
    Return a dictionary of data which is a list of tuples.
    '''
    d_of_data = {}

    d_of_data['score'] = get_score_data(pose, structure_name)
    d_of_data['residue_average_score'] = get_residue_average_score_data(pose, structure_name)
    d_of_data['residue_score'] = get_residue_score_data(pose, structure_name)
    
    d_of_data['num_buried_unsatisfied_hbonds'] = get_num_buried_unsatisfied_hbonds_data(pose, structure_name)
    d_of_data['average_num_buried_unsatisfied_hbonds'] = get_average_num_buried_unsatisfied_hbonds_data(pose, structure_name)
    
    d_of_data['num_over_saturated_hbond_acceptors'] = get_num_over_saturated_hbond_acceptors_data(pose, structure_name)
    d_of_data['average_num_over_saturated_hbond_acceptors'] = get_average_num_over_saturated_hbond_acceptors_data(pose, structure_name)

    d_of_data['hydrophobic_sasa'] = get_hydrophobic_sasa_data(pose, structure_name)
    d_of_data['average_hydrophobic_sasa'] = get_average_hydrophobic_sasa_data(pose, structure_name)
   
    d_of_data['holes_score'] = get_holes_score_data(pose, structure_name)

    return d_of_data

def get_metrics(data_path, input_path, num_jobs, job_id):
    '''Get rosetta metrics for all structures in the
    input path.
    '''
    pyrosetta.init('-ignore_unrecognized_res true')

    for i, f in enumerate(os.listdir(input_path)):
        if i % num_jobs != job_id:
            continue

        if f.endswith('.pdb') or f.endswith('.pdb.gz'):
            pose = rosetta.core.pose.Pose()
            rosetta.core.import_pose.pose_from_file(pose, os.path.join(input_path, f))

            d_of_data = get_metrics_for_one_pose(pose, f.split('.')[0])
            
            for key in d_of_data.keys():
                append_to_csv(d_of_data[key], os.path.join(data_path, '{0}.csv'.format(key)))


if __name__ == '__main__':
    
    data_path = sys.argv[1]
    input_path = sys.argv[2]
    
    num_jobs = 1
    job_id = 0
    if len(sys.argv) > 4:
      num_jobs = int(sys.argv[3])
      job_id = int(sys.argv[4]) - 1
    
    get_metrics(data_path, input_path, num_jobs, job_id)
