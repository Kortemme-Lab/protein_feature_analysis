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


class StructuralHomologFeature(Feature):
  '''Analyze features of structrual homologs.'''
  
  def __init__(self, dssp_path):
    super().__init__()
    self.superfamilies = []
    self.dssp_path = dssp_path

  def extract(self, input_path, output_path, total_num_threads=1, my_id=0):
    '''Extract structrual homolog features from structures in the input path.
    The input data should be stored in .pml files of superposed homologous 
    superfamilies from the CATH database.
    '''

    for f in self.list_my_jobs(input_path, total_num_threads, my_id):
      if f.endswith('.pml'):
      
        # Make a scratch directory

        scratch_path = os.path.join(output_path, f[0:-4])
        if not os.path.exists(scratch_path):
          os.mkdir(scratch_path)

        # Extract features from one file  

        self.extract_from_one_file(os.path.join(input_path, f), scratch_path)

  def extract_from_one_file(self, pml_file, scratch_path):
    '''Extract structrual homolog features from a .pml file of superposed
    homologous superfamilies from the CATH database.
    '''
    self.superfamilies.append([])
    candidate_proteins = []

    with open(pml_file, 'r') as f:
      while True:
        line = f.readline()
        if not line: break
       
        # Read one structure

        if line.strip().startswith('cmd.read_pdbstr'):
          pdb_lines = [line.strip()[19:].strip('\\')]
          pdb_id = ''

          while True:  
            line = f.readline()
            if line.strip().startswith('"""'):
              pdb_id = line.strip()[5:12]
              break

            pdb_line = line.strip().strip('\\')
            if len(pdb_line) > 17:
              pdb_line = pdb_line[0:16] + ' ' + pdb_line[17:] # Remove all altLoc flags
            
            pdb_lines.append(pdb_line) # Remove all altLoc flags
         
          # Make a pdb file of the structure for DSSP analysis
         
          structure = self.structure_from_pdb_string('\n'.join(pdb_lines), pdb_id)

          # Store structures without chain breaks

          if len(topology.find_structure_chain_breaks(structure)) == 0:
            structure_path = os.path.join(scratch_path, pdb_id + '.pdb')

            io = PDB.PDBIO()
            io.set_structure(structure)
            io.save(structure_path)

            candidate_proteins.append({'structure' : structure, 'path' : structure_path})

    for p in candidate_proteins:
      try:  
        self.find_secondary_structures(p)
      except:
        continue
      self.superfamilies[-1].append(p) # Add a protein to a superfamily if there's no exception

  def find_secondary_structures(self, protein_dict):
    '''Find secondary structures of a protein.
    
    Arguements:
     - protein_dict - a dictionary to store informations of a protein
    '''
    protein_dict['dssp_dict'], protein_dict['dssp_key_map'] = \
      secondary_structures.make_dssp_dict(protein_dict['path'], self.dssp_path)

    protein_dict['ss_list'], protein_dict['sheet_list'] = \
      secondary_structures.pack_dssp_dict_into_ss_list(protein_dict['structure'][0],
            protein_dict['dssp_dict'], protein_dict['dssp_key_map'])


  ############## Functions for calculating, saving, loading and visualizing features #################

  def calc_ss_composition_features(self):
    '''Calculate the secondary structure compositions.'''
    self.ss_composition_features = []

    for sf in self.superfamilies:
      for p in sf:
        d = {}
        d['num_alpha_helix'] = len([ss for ss in p['ss_list'] if isinstance(ss, secondary_structures.AlphaHelix)]) 
        d['num_strand'] = len([ss for ss in p['ss_list'] if isinstance(ss, secondary_structures.BetaStrand)]) 
        d['num_3-10_helix'] = len([ss for ss in p['ss_list'] if isinstance(ss, secondary_structures.Loop) and ss.type == 'G']) 
        d['num_PI_helix'] = len([ss for ss in p['ss_list'] if isinstance(ss, secondary_structures.Loop) and ss.type == 'I']) 
        d['num_turn'] = len([ss for ss in p['ss_list'] if isinstance(ss, secondary_structures.Loop) and ss.type == 'T'])
        d['num_sheet'] = len(p['sheet_list'])
        d['pdb_id'] = p['structure'].get_id()
 
        self.ss_composition_features.append(d)

  def save_ss_composition_features(self, data_path):
    '''Save the secondary structure compositions.'''
    self.calc_ss_composition_features()
    df = pd.DataFrame(data=self.ss_composition_features)
    
    self.append_to_csv(df, os.path.join(data_path, 'ss_composition_features.csv'))

  def load_ss_composition_features(self, data_path):
    '''Load the secondary structure compositions.'''
    df = pd.read_csv(os.path.join(data_path, 'ss_composition_features.csv'), header=None)
   
    self.ss_composition_features = []
    for index, row in df.iterrows():
      self.ss_composition_features.append({'num_3-10_helix':row[0], 'num_PI_helix':row[1],
          'num_alpha_helix':row[2], 'num_sheet':row[3], 'num_strand':row[4], 'num_turn':row[5], 'pdb_id':row[6]})
 
  def visualize_ss_composition_features(self, vis_type='num_alpha_helix'):
    '''Visualize the secondary structure compositions.'''
    data = [d[vis_type] for d in self.ss_composition_features]
    hist, bin_edges = np.histogram(data, bins=np.arange(max(data)))

    plt.bar(bin_edges[0:-1], hist, width=1, edgecolor='black')
    plt.show()

  def calc_ss_sizes_features(self):
    '''Calculate sizes of secondary structures.'''
    self.ss_sizes_features = []

    for sf in self.superfamilies:
      for p in sf:
         
          # Linear structures

          for i in range(len(p['ss_list'])):

            if isinstance(p['ss_list'][i], secondary_structures.AlphaHelix):
              self.ss_sizes_features.append({'ss':'alpha_helix', 'type':'-',
                  'size' : len(p['ss_list'][i].residue_list)})

            elif isinstance(p['ss_list'][i], secondary_structures.BetaStrand):
              self.ss_sizes_features.append({'ss':'beta_strand', 'type':'-',
                  'size' : len(p['ss_list'][i].residue_list)})

            elif isinstance(p['ss_list'][i], secondary_structures.Loop):
              l_len = len(p['ss_list'][i].residue_list)
               
              j = i
              while j + 1 < len(p['ss_list']) and isinstance(p['ss_list'][j + 1], secondary_structures.Loop):
                j += 1
                l_len += len(p['ss_list'][j].residue_list)

              l_type = '-'
              if i == 0 or j == len(p['ss_list']) - 1:
                pass
              elif isinstance(p['ss_list'][i - 1], secondary_structures.AlphaHelix) \
                   and isinstance(p['ss_list'][j + 1], secondary_structures.AlphaHelix):
                l_type = 'alpha-alpha'
              elif isinstance(p['ss_list'][i - 1], secondary_structures.AlphaHelix) \
                   and isinstance(p['ss_list'][j + 1], secondary_structures.BetaStrand):
                l_type = 'alpha-beta'
              elif isinstance(p['ss_list'][i - 1], secondary_structures.BetaStrand) \
                   and isinstance(p['ss_list'][j + 1], secondary_structures.AlphaHelix):
                l_type = 'beta-alpha'
              elif isinstance(p['ss_list'][i - 1], secondary_structures.BetaStrand) \
                   and isinstance(p['ss_list'][j + 1], secondary_structures.BetaStrand):
                l_type = 'beta-beta'

              self.ss_sizes_features.append({'ss':'loop', 'type':l_type,
                  'size' : l_len})

              i = j + 1

          # Beta sheets

          for sheet in p['sheet_list']:

            self.ss_sizes_features.append({'ss':'beta_sheet', 'type':sheet.type,
                'size':len(sheet.strand_list)})

  def save_ss_sizes_features(self, data_path):
    '''Save sizes of secondary structures.'''
    self.calc_ss_sizes_features()
    df = pd.DataFrame(data=self.ss_sizes_features)
    
    self.append_to_csv(df, os.path.join(data_path, 'ss_sizes_features.csv'))

  def load_ss_sizes_features(self, data_path):
    '''Load sizes of secondary structures.'''
    df = pd.read_csv(os.path.join(data_path, 'ss_sizes_features.csv'), header=None)
   
    self.ss_sizes_features = []
    for index, row in df.iterrows():
      self.ss_sizes_features.append({'size':row[0], 'ss':row[1],
          'type':row[2]})
 
  def visualize_ss_sizes_features(self, ss='alpha_helix', s_type='-', fig_save_path=None):
    '''Visualize sizes of secondary structures.'''
    data = [d['size'] for d in self.ss_sizes_features if d['ss'] == ss and d['type'] == s_type]
    hist, bin_edges = np.histogram(data, bins=np.arange(max(data)))

    plt.bar(bin_edges[0:-1], hist, width=1, edgecolor='black')
    plt.ylabel('Number of structures')
    plt.xlabel('Length of ' + ss + ' - ' + s_type)
    
    if fig_save_path:
      plt.savefig(os.path.join(fig_save_path, ss + '-' + s_type + '.png'))
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
            self.beta_sheet_features.append({
                'edge_cur1' : d['edge_cur1'],
                'edge_cur2' : d['edge_cur2'],
                'edge_cur3' : d['edge_cur3'],
                'edge_cur4' : d['edge_cur4'],
                'gauss_cur' : d['gauss_cur'],
                'mismatch' : d['mismatch'], 'path' : p['path'], 
                'phi' : np.degrees(geometry.get_phi(res.get_parent(), res)), 
                'psi' : np.degrees(geometry.get_psi(res.get_parent(), res)), 
                'side' : d['side'], 'terminal' : d['terminal'], 'type' : sheet.type})

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
          'edge_cur1':row[0], 'edge_cur2':row[1], 'edge_cur3':row[2],'edge_cur4':row[3],
          'gauss_cur':row[4], 'mismatch':row[5], 'path':row[6],
          'phi':row[7], 'psi':row[8], 'side':row[9], 'terminal':row[10], 'type':row[11]})

  def visualize_beta_sheet_features(self, sheet_type, feature1, feature2='', 
          side_threshold=(0, 10000), terminal_threshold=(0, 10000), 
          mismatch_threshold=(0, 10000), fig_save_path=None):
    '''Visualize beta sheet features.'''

    def get_data(feature):
      data = [] 
      for d in self.beta_sheet_features: 
        if d['type'] == sheet_type \
            and side_threshold[0] <= d['side'] <= side_threshold[1] \
            and terminal_threshold[0] <= d['terminal'] <= terminal_threshold[1] \
            and mismatch_threshold[0] <= d['mismatch'] <= mismatch_threshold[1] :
          
          keys = [k for k in d.keys() if k.startswith(feature)]
          for k in keys:
            if d[k] != 'null':
              data.append(float(d[k]))
      
      return data

    data1 = get_data(feature1)

    # Visualize one feature
    
    if feature2 == '':

      # Calculate mean and standard deviation
      
      print("mean = {0}, std = {1}, num_data = {2}".format(
          np.mean(data1), np.std(data1), len(data1)))
      
      # Make histograms
      
      step = (max(data1) - min(data1)) / 100
      hist, bin_edges = np.histogram(data1, bins=np.arange(min(data1), max(data1), step))

      plt.clf()
      plt.bar(bin_edges[0:-1], hist, width=step, edgecolor='black')
      plt.ylabel('Number of residues')
      plt.xlabel('-'.join([sheet_type, feature1]))
  
    # Visualize two features

    else:
      
      data2 = get_data(feature2)

      # Draw a heat map 
    
      heatmap, xedges, yedges = np.histogram2d(data1, data2, bins=(128,128))
      extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
      plt.imshow(np.transpose(heatmap), extent=extent, origin='lower', aspect='auto')

    # Save or plot

    if fig_save_path:
      plt.savefig(os.path.join(fig_save_path, '-'.join(['beta_sheet', sheet_type, feature1, feature2]) + '.png'))
    else:
      plt.show()
