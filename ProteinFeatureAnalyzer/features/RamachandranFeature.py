import os

import matplotlib
matplotlib.use('TkAgg') 
import numpy as np
import matplotlib.pyplot as plt
import Bio.PDB as PDB

from .Feature import Feature


class RamachandranFeature(Feature):
  '''Analyze the phi/psi torsion distributions of proteins.'''
  
  def __init__(self):
    super().__init__()

  def extract(self, input_path):
    '''Extract phi, psi angles from structures in the input path.'''
    for f in os.listdir(input_path):
      if f.endswith('.pdb'):
        self.extract_from_one_file(os.path.join(input_path, f))

  def extract_from_one_file(self, pdb_file):
    '''Extract phi, psi angles from a pdb_file.'''
    structure = self.structure_from_pdb_file(pdb_file)

    for model in structure:
      for chain in model:
        for residue in chain:
          try:
            feature_dict = {'phi' : self.get_phi(chain, residue),
                            'psi' : self.get_psi(chain, residue)}
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
