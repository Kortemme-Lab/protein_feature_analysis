import os
import json

import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import pandas as pd
import Bio.PDB as PDB

from .Feature import Feature
from . import topology
from . import geometry
from . import secondary_structures
from . import data_loading


class ParametricDesignFeature(Feature):
  '''Analyze features for parametric design.'''
  
  def __init__(self, dssp_path):
    super().__init__()
    self.superfamilies = []
    self.dssp_path = dssp_path

  def extract(self, input_path, output_path, total_num_threads=1, my_id=0):
    '''Extract structrual homolog features from structures in the input path.
    The input data should be stored in .pml files of superposed homologous 
    superfamilies from the CATH database.
    '''

    self.superfamilies = data_loading.load_data_from_cath_pmls(input_path, output_path,
                            self.list_my_jobs(input_path, total_num_threads, my_id), self.dssp_path)

  ############## Functions for calculating, saving, loading and visualizing features #################

  def calc_alpha_helix_features(self):
    '''Calculate features of alpha helices.'''
    self.alpha_helix_features = []

    for sf in self.superfamilies:
      for p in sf:
        for helix in p['ss_list']:
          if not isinstance(helix, secondary_structures.AlphaHelix):
              continue

          self.alpha_helix_features.append({'angles': helix.angles, 'dihedrals': helix.dihedrals})

  def save_alpha_helix_features(self, data_path):
    '''Save alpha helix features.'''
    self.calc_alpha_helix_features()
    df = pd.DataFrame(data=self.alpha_helix_features)
    
    self.append_to_csv(df, os.path.join(data_path, 'alpha_helix_features.csv'))

  def load_alpha_helix_features(self, data_path):
    '''Load alpha helix features.'''
    df = pd.read_csv(os.path.join(data_path, 'alpha_helix_features.csv'), header=None)
   
    self.alpha_helix_features = []
    for index, row in df.iterrows():
        self.alpha_helix_features.append({'angles':json.loads(row[0]), 'dihedrals':json.loads(row[1])})

  def visualize_alpha_helix_features(self, feature, fig_save_path=None):
    '''Visualized the alpha helix features.'''
    
    def get_data(feature):
      data = []

      for d in self.alpha_helix_features:
        data += d[feature]

      return [np.degrees(x) for x in data]
    
    data = get_data(feature)

    # Calculate mean and standard deviation
    
    print("mean = {0}, std = {1}, num_data = {2}".format(
        np.mean(data), np.std(data), len(data)))
    
    # Make histograms
    
    step = (max(data) - min(data)) / 100
    hist, bin_edges = np.histogram(data, bins=np.arange(min(data), max(data), step))

    plt.clf()
    plt.bar(bin_edges[0:-1], hist, width=step, edgecolor='black')
    plt.ylabel('Number of data')
    plt.xlabel(feature)
    
    # Save or plot

    if fig_save_path:
      plt.savefig(os.path.join(fig_save_path, '-'.join(['alpha_helix', feature]) + '.svg'))
    else:
      plt.show()
  
  def calc_beta_sheet_features(self):
    '''Calculate features of beta sheets.'''
    self.beta_sheet_features = []

    for sf in self.superfamilies:
      for p in sf:
        for sheet in p['sheet_list']:
          for res in sheet.graph.nodes():

            d = sheet.graph.node[res]
            if d['cylinder_radius'] != 'null':
                self.beta_sheet_features.append({'cylinder_radius' : d['cylinder_radius'],
                    'cylinder_strand_angle' : d['cylinder_strand_angle'],
                    'cylinder_fitting_rmsd' : d['cylinder_fitting_rmsd'],
                    'folding_angle' : d['folding_angle']})

  def save_beta_sheet_features(self, data_path):
    '''Save beta sheet features.'''
    self.calc_beta_sheet_features()
    df = pd.DataFrame(data=self.beta_sheet_features)
    
    self.append_to_csv(df, os.path.join(data_path, 'beta_sheet_features.csv'))

  def load_beta_sheet_features(self, data_path):
    '''Load beta sheet features.'''
    df = pd.read_csv(os.path.join(data_path, 'beta_sheet_features.csv'), header=None)
   
    self.beta_sheet_features = []
    for index, row in df.iterrows():
        self.beta_sheet_features.append({
            'cylinder_fitting_rmsd' : row[0],
            'cylinder_radius' : row[1],
            'cylinder_strand_angle' : row[2],
            'folding_angle' : row[3]})

  def visualize_beta_sheet_features(self, feature, fig_save_path=None):
    '''Visualized the beta sheet features features.'''
    
    def get_data(feature):
      if feature == 'cylinder_curvature':
          return [1 / d['cylinder_radius'] for d in self.beta_sheet_features]
      elif feature in ['cylinder_strand_angle', 'folding_angle']:
          return [np.degrees(d[feature]) for d in self.beta_sheet_features]
      else:
          return [d[feature] for d in self.beta_sheet_features]
    
    data = get_data(feature)

    # Calculate mean and standard deviation
    
    print("mean = {0}, std = {1}, num_data = {2}".format(
        np.mean(data), np.std(data), len(data)))
    
    # Make histograms
    
    step = (max(data) - min(data)) / 100
    hist, bin_edges = np.histogram(data, bins=np.arange(min(data), max(data), step))

    plt.clf()
    plt.bar(bin_edges[0:-1], hist, width=step, edgecolor='black')
    plt.ylabel('Number of data')
    plt.xlabel(feature)
    
    # Save or plot

    if fig_save_path:
      plt.savefig(os.path.join(fig_save_path, '-'.join(['beta_sheet', feature]) + '.svg'))
    else:
      plt.show()
  

