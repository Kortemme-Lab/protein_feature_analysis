import os

import matplotlib
matplotlib.use('TkAgg') 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Bio.PDB as PDB

from .Feature import Feature
from . import Geometry


class RamachandranFeature(Feature):
  '''Analyze the phi/psi torsion distributions of proteins.'''
  
  def __init__(self):
    super().__init__()

  def extract(self, input_path, total_num_threads=1, my_id=0):
    '''Extract phi, psi angles from structures in the input path.'''
    for f in self.list_my_jobs(input_path, total_num_threads, my_id):
      if f.endswith('.pdb'):
        self.extract_from_one_file(os.path.join(input_path, f))

  def extract_from_one_file(self, pdb_file):
    '''Extract phi, psi angles from a pdb_file.'''
    structure = self.structure_from_pdb_file(pdb_file)

    for model in structure:
      for chain in model:
        for residue in chain:
          try:
            feature_dict = {'phi' : Geometry.get_phi(chain, residue),
                            'psi' : Geometry.get_psi(chain, residue)}
            self.feature_list.append(feature_dict)
          except:
            pass

  def visualize(self):
    '''Visualize the feature statistics.'''
    phis = [ d['phi'] for d in self.feature_list ]
    psis = [ d['psi'] for d in self.feature_list ]

    plt.scatter(phis, psis)

    plt.axis([- np.pi, np.pi, - np.pi, np.pi])
    plt.show()

  def save(self, data_path):
    '''Save the data into a csv file.'''
    data = [ (d['phi'], d['psi']) for d in self.feature_list ]
    df = pd.DataFrame(data=data, columns=['phi', 'psi'])
    
    self.append_to_csv(df, os.path.join(data_path, 'rama_features.csv'))

  def load(self, data_path):
    '''Load data from a csv file.'''
    df = pd.read_csv(os.path.join(data_path, 'rama_features.csv'), header=None)
    
    for index, row in df.iterrows():
      self.feature_list.append({'phi':row[0], 'psi':row[1]})
