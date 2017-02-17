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

