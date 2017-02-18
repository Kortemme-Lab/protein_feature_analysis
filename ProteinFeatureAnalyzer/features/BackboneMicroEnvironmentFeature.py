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

        # Get the relative orientation of the second residue

        s_matrix2, origin2 = Geometry.get_residue_stub_matrix(res2)
        rot_matrix = np.dot(s_matrix.T, s_matrix2) # Rotation matrix in the frame of the first residue
        feature_dict['theta_x'], feature_dict['theta_y'], feature_dict['theta_z'] = \
                Geometry.rotation_matrix_to_euler_angles(np.array(rot_matrix))

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

    # Postions of N and C

    n_ca_c_angle = 110.86 * np.pi / 180

    fig =plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(X, Y, Z, c='green')
    ax.quiver([0], [0], [0], [1.32869], [0], [0], color='blue')
    ax.quiver([0], [0], [0], [1.52326 * np.cos(n_ca_c_angle)], [1.52326 * np.sin(n_ca_c_angle)], [0], color='red')

    plt.show()

  def plot_shift_length_histogram(self):
    '''Plot a histogram of the lengths of translational shifts.'''

    lengths = [np.linalg.norm(d['shift']) for d in self.feature_list]
    hist, bin_edges = np.histogram(lengths, bins=0.5 * np.arange(20))

    plt.bar(bin_edges[0:-1] - 0.25, hist, width=0.5, edgecolor='black')
    plt.show()

  def scatter_plot_two_features(self, feature1_l, feature2_l, axis=None):
    '''Make a scatter plot of two features. feature1_l and feature2_l
    are lambda expressions for pick a feature.
    '''

    f1 = [feature1_l(d) for d in self.feature_list]
    f2 = [feature2_l(d) for d in self.feature_list]

    plt.scatter(f1, f2)
    if axis:
      plt.axis(axis)
    plt.show()
