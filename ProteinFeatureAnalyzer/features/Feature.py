import os

import numpy as np
import Bio.PDB as PDB


class Feature:
  '''Base class for features.'''

  def __init__(self):
    self.feature_list = []

  def list_my_jobs(self, input_path, total_num_threads, my_id):
    '''List all the inputs that should be hanled by the running thread.'''
    all_jobs = os.listdir(input_path)
    return [all_jobs[i] for i in range(len(all_jobs)) if i % total_num_threads == my_id]


  def structure_from_pdb_file(self, file_path):
    '''Read the structure stored in a PDB file.'''
    parser = PDB.PDBParser()
    return parser.get_structure('', file_path)

  def get_phi(self, chain, residue):
    '''Calculate the phi torsion of a residue.'''
    
    # Get the previous residue

    res_id = residue.get_id()
    prev_res = chain[res_id[1] - 1]
    
    prev_flag = prev_res.get_id()[0]
    if prev_flag == 'W' or prev_flag.startswith('H_'):
      raise Exception('Hetero residue type!')
   
    # Calculate the torsion

    c_prev = prev_res['C'].get_vector()
    n = residue['N'].get_vector()
    ca = residue['CA'].get_vector()
    c = residue['C'].get_vector()

    return PDB.calc_dihedral(c_prev, n, ca, c) 
  
  def get_psi(self, chain, residue):
    '''Calculate the psi torsion of a residue.'''
    
    # Get the next residue

    res_id = residue.get_id()
    next_res = chain[res_id[1] + 1]
    
    next_flag = next_res.get_id()[0]
    if next_flag == 'W' or next_flag.startswith('H_'):
      raise Exception('Hetero residue type!')
    
    # Calculate the torsion

    n = residue['N'].get_vector()
    ca = residue['CA'].get_vector()
    c = residue['C'].get_vector()
    n_next = next_res['N'].get_vector()

    return PDB.calc_dihedral(n, ca, c, n_next) 
