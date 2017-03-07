import os

import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
import pandas as pd
import Bio.PDB as PDB

from .Feature import Feature
from . import topology
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

