import os

import numpy as np
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
  
    for model in structure:
      nearest_nb_list = Geometry.get_nearest_nonbonded_residues(model)
   
      for res1, res2 in nearest_nb_list:
        feature_dict ={}
        
        # Get the torsions
        
        try:
          feature_dict['phi1'] = Geometry.get_phi(res1.get_parent(), res1) 
          feature_dict['psi1'] = Geometry.get_psi(res1.get_parent(), res1) 
          feature_dict['phi2'] = Geometry.get_phi(res2.get_parent(), res1) 
          feature_dict['psi2'] = Geometry.get_psi(res2.get_parent(), res1) 
        except:
          continue

        # Get the relative position of the second residue 

        s_matrix, origin = Geometry.get_residue_stub_matrix(res1)
        shift_global = res2['CA'].get_coord() - res1['CA'].get_coord()
        feature_dict['shift'] = np.matmul(np.array(s_matrix.T), np.array(shift_global))

        self.feature_list.append(feature_dict)

    #print(self.feature_list)

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

  def plot_shifts(self):
    '''Plot the distribution of the translational shifts from the 
    CA atom of the first residue to the CA atom of the second residue.
    '''

    # Data points

    X = [d['shift'][0] for d in self.feature_list]
    Y = [d['shift'][1] for d in self.feature_list]
    Z = [d['shift'][2] for d in self.feature_list]

    # Postions of N and ca

    fig =plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(X, Y, Z, c='green')
    plt.show()
