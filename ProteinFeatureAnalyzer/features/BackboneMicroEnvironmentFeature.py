import os

import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

from .Feature import Feature
from . import Geometry


class BackboneMicroEnvironmentFeature(Feature):
  '''The BackboneMicroEnvironmentFeature analyzes the micro environments
  of backbones of each residue. The micro environemt is formed by the residue
  it self and its nearest non-connected residue.'''
  
  def __init__(self):
    super().__init__()
  
  def extract(self, input_path, total_num_threads=1, my_id=0):
    '''Extract features from structures in the input path.'''
    for f in self.list_my_jobs(input_path, total_num_threads, my_id):
      if f.endswith('.pdb'):
        self.extract_from_one_file(os.path.join(input_path, f))

  def extract_from_one_file(self, pdb_file):
    structure = self.structure_from_pdb_file(pdb_file)
   
    nearest_nb_list = Geometry.get_nearest_nonbonded_residues(structure)
    self.plot_nearst_nonbonded_list(nearest_nb_list)

  def visualize(self):
    pass

  def plot_nearst_nonbonded_list(self, nearest_nb_list):
    '''A debugging function to print the nearest nonbonded residue list.'''
    X = [pair[0]['CA'].get_coord()[0] for pair in nearest_nb_list]
    Y = [pair[0]['CA'].get_coord()[1] for pair in nearest_nb_list]
    Z = [pair[0]['CA'].get_coord()[2] for pair in nearest_nb_list]
    
    U = [pair[1]['CA'].get_coord()[0] - pair[0]['CA'].get_coord()[0] for pair in nearest_nb_list]
    V = [pair[1]['CA'].get_coord()[1] - pair[0]['CA'].get_coord()[1] for pair in nearest_nb_list]
    W = [pair[1]['CA'].get_coord()[2] - pair[0]['CA'].get_coord()[2] for pair in nearest_nb_list]

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.quiver(X, Y, Z, U, V, W)
    plt.show()
