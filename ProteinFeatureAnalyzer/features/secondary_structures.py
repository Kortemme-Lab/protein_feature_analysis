import subprocess

import numpy as np
import Bio.PDB as PDB

from . import geometry


def make_dssp_dict(pdb_path, dssp_path):
  '''Make a dictionary that stores the DSSP information of a
  protein. The keys are (chain_id, res_id).

  Return values:
    - a DSSP information dictionary
    - a map from dssp sequential number to keys
  '''
  # Get the DSSP outputs

  p = subprocess.Popen([dssp_path, pdb_path], universal_newlines=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  out, err = p.communicate()
  
  # Report error
  
  if err.strip():
    warnings.warn(err)
    if not out.strip():
      raise Exception('DSSP failed to produce an output')

  # Parse the DSSP output

  lines = out.split('\n')
  i = 0

  # Ignore the headers

  while len(lines[i].split()) < 2 or lines[i].split()[1] != 'RESIDUE':
    i += 1

  # Parse line by line

  dssp_dict = {}
  key_map = {}

  for l in lines[i + 1:]:
    if len(l.strip()) == 0: continue
    if l[9] == " ":
            # Skip -- missing residue
            continue

    seq_num = int(l[0:5])   # Sequential residue number,including chain breaks as extra residues
    res_id = int(l[5:10])   # Residue ID in PDB file
    chain_id = l[11]        # Chain ID
    ss = l[16]              # Secondary structure
    bbl1 = l[23]            # Beta bridge label 1
    bbl2 = l[24]            # Beta bridge label 2
    bp1 = int(l[26:29])     # Beta bridge partner 1
    bp2 = int(l[30:33])     # Beta bridge partner 2
    bsl = l[33]             # Beta sheet label

    key_map[seq_num] = (chain_id, res_id)
    dssp_dict[(chain_id, res_id)] = {'seq_num':seq_num, 'ss':ss,
            'bbl1':bbl1, 'bbl2':bbl2, 'bp1':bp1, 'bp2':bp2, 'bsl':bsl}

  return dssp_dict, key_map


def pack_dssp_dict_into_ss_list(model, dssp_dict, key_map):
  '''Pack a dssp dictionary into a list of secondary structures
  and a list of beta sheets.
  '''
 
  # map to translate secondary structures

  ss_map = {'H':'H', 'B':'B', 'E':'E', 'G':'G', 'I':'I', 'T':'T'}

  def same_ss(dssp_dict, key1, key2):
    '''Return if two residues have save secondary structures.
    Bend are regarded as normal loops.
    Note that for beta strands, they should be in the same sheet
    to be considered the same.
    '''
    return ss_map.setdefault(dssp_dict[key1]['ss'], ' ') == ss_map.setdefault(dssp_dict[key2]['ss'], ' ') \
            and dssp_dict[key1]['bsl'] == dssp_dict[key2]['bsl']


  # Get the list of secondary structures
  
  ss_list = []

  for chain in model:
    c_id = chain.get_id()

    residue_list = [r for r in chain]

    start = 0
    while start < len(residue_list):
      start_r_id = residue_list[start].get_id()[1]
      start_ss = dssp_dict[(c_id, start_r_id)]['ss']

      # Find the stop position of the secondary structure

      stop = start + 1
      while stop < len(residue_list):
        stop_r_id = residue_list[stop].get_id()[1]
        
        if not same_ss(dssp_dict, (c_id, start_r_id), (c_id, stop_r_id)):
          break
        stop += 1

      stop -= 1

      # Add the secondary structure into a list

      ss_residue_list = [residue_list[i] for i in range(start, stop + 1)] 

      if 'H' == start_ss:
        ss_list.append(AlphaHelix(ss_residue_list))
      
      elif 'E' == start_ss or 'B' == start_ss:
        ss_list.append(BetaStrand(ss_residue_list, dssp_dict[(c_id, start_r_id)]['bsl']))
      
      else:
        ss_list.append(Loop(ss_residue_list, ss_map.setdefault(start_ss, ' ')))

      start = stop + 1

  # Get the list of beta sheets

  sheet_list = []
  sheet_id_set = set()
  for ss in ss_list:
    if isinstance(ss, BetaStrand):
      sheet_id_set.add(ss.sheet_id)

  for sheet_id in sheet_id_set:  
    sheet_list.append(BetaSheet(sheet_id, 
        [ss for ss in ss_list if isinstance(ss, BetaStrand) and ss.sheet_id == sheet_id],
        dssp_dict))

  return ss_list, sheet_list


class SecondaryStructure:
  '''Base class for secondary structures.'''
  def __init__(self):
    pass


class AlphaHelix(SecondaryStructure):
  '''Class that represents an alpha helix.'''
  def __init__(self, residue_list):
    self.residue_list = residue_list

  def parameterize(self):
    '''Get the Ca representation of the helix and get the bond lengths,
    bond angles and torsions of the Ca representation. There will be n-1 bonds,
    n-2 angles and n-3 torsions.
    '''
    cas = [r['CA'] for r in self.residue_list]
    cas_v = [a.get_vector() for a in cas]
    
    self.bond_lengths = []
    self.bond_angles = []
    self.bond_torsions = []
    
    # Get bond lengths
    
    for i in range(len(cas) - 1):
      self.bond_lengths.append(cas[i + 1] - cas[i])

    # Get bond angles

    for i in range(1, len(cas) - 1):
      self.bond_angles.append(np.degrees(PDB.calc_angle(cas_v[i - 1],
          cas_v[i], cas_v[i + 1])))

    # Get torsions

    for i in range(1, len(cas) - 2):
      self.bond_torsions.append(np.degrees(PDB.calc_dihedral(cas_v[i - 1],
          cas_v[i], cas_v[i + 1], cas_v[i + 2])))

class BetaStrand(SecondaryStructure):
  '''Class that represents a beta strand.'''
  def __init__(self, residue_list, sheet_id):
    self.residue_list = residue_list
    self.sheet_id = sheet_id

  def calc_bb_threading_features(self):
    '''Calculate backbone threading features.'''
    self.ca_n_internal_coord_n_term = []
    self.ca_c_internal_coord_n_term = []
    self.ca_n_internal_coord_mid = []
    self.ca_c_internal_coord_mid = []
    self.ca_n_internal_coord_c_term = []
    self.ca_c_internal_coord_c_term = []

    if len(self.residue_list) < 3:
        return

    # Calculate the N terminal features

    self.ca_n_internal_coord_n_term.append(
            geometry.cartesian_coord_to_internal_coord(
            self.residue_list[2]['CA'].get_coord(),
            self.residue_list[1]['CA'].get_coord(),
            self.residue_list[0]['CA'].get_coord(),
            self.residue_list[0]['N'].get_coord()))

    self.ca_c_internal_coord_n_term.append( 
            geometry.cartesian_coord_to_internal_coord(
            self.residue_list[2]['CA'].get_coord(),
            self.residue_list[1]['CA'].get_coord(),
            self.residue_list[0]['CA'].get_coord(),
            self.residue_list[0]['C'].get_coord()))

    # Calculate the middle features

    for i in range(1, len(self.residue_list) - 1):

      self.ca_n_internal_coord_mid.append(
              geometry.cartesian_coord_to_internal_coord(
              self.residue_list[i - 1]['CA'].get_coord(),
              self.residue_list[i + 1]['CA'].get_coord(),
              self.residue_list[i]['CA'].get_coord(),
              self.residue_list[i]['N'].get_coord()))

      self.ca_c_internal_coord_mid.append(
              geometry.cartesian_coord_to_internal_coord(
              self.residue_list[i + 1]['CA'].get_coord(),
              self.residue_list[i - 1]['CA'].get_coord(),
              self.residue_list[i]['CA'].get_coord(),
              self.residue_list[i]['C'].get_coord()))

    # Calculate the C terminal features

    self.ca_n_internal_coord_c_term.append(
            geometry.cartesian_coord_to_internal_coord(
            self.residue_list[-3]['CA'].get_coord(),
            self.residue_list[-2]['CA'].get_coord(),
            self.residue_list[-1]['CA'].get_coord(),
            self.residue_list[-1]['N'].get_coord()))

    self.ca_c_internal_coord_c_term.append(
            geometry.cartesian_coord_to_internal_coord(
            self.residue_list[-3]['CA'].get_coord(),
            self.residue_list[-2]['CA'].get_coord(),
            self.residue_list[-1]['CA'].get_coord(),
            self.residue_list[-1]['C'].get_coord()))

class BetaSheet(SecondaryStructure):
  '''Class that represents a beta sheet.'''
  def __init__(self, sheet_id, strand_list, dssp_dict):
    self.sheet_id = sheet_id
    self.strand_list = strand_list
    self.init_sheet_network(dssp_dict)
    self.determine_sheet_type()

  def init_sheet_network(self, dssp_dict):
    '''Initialize a network that represents the beta sheet.'''
    self.sheet_network = []

    # Initialize a none-connected network

    for strand in self.strand_list:
      for r in strand.residue_list:
        self.sheet_network.append({'residue':r, 'prev':None, 'next':None, 'bps':[]})

    # Initialize connections within strands

    for i in range(len(self.sheet_network) - 1):
      node1 = self.sheet_network[i]
      node2 = self.sheet_network[i + 1]

      if node1['residue'].get_parent().get_id() != node2['residue'].get_parent().get_id():
        continue

      if node1['residue'].get_id()[1] + 1 == node2['residue'].get_id()[1]:
        node1['next'] = i + 1
        node2['prev'] = i

    # Initialize a map from the DSSP sequential numbers to indices of the sheet_network

    id_map = {}

    for i, node in enumerate(self.sheet_network):
      seq_num = dssp_dict[(node['residue'].get_parent().get_id(),
                           node['residue'].get_id()[1])]['seq_num']
      id_map[seq_num] = i

    # Initialize connections between strands
   
    for node in self.sheet_network:
      res_id = (node['residue'].get_parent().get_id(), node['residue'].get_id()[1])
      bp1 = dssp_dict[res_id]['bp1']
      bp2 = dssp_dict[res_id]['bp2']

      if bp1 > 0:
        node['bps'].append(id_map[bp1])
      if bp2 > 0:
        node['bps'].append(id_map[bp2])

  def determine_sheet_type(self):
    '''Determine if the sheet is parallel, antiparallel or mixed.'''
    num_parallel = 0
    num_antiparallel = 0

    i = 0
    while i < len(self.sheet_network):
      
      # Find the start of this strand

      strand_start = i

      while self.sheet_network[i]['next']:
        i += 1
      i += 1

      # If the strand only has one residue, take add a antiparallel strand

      if self.sheet_network[strand_start]['next'] == None:
        num_antiparallel += 1
        continue

      # Determine if the partner strands are parallel to this strand

      for j in self.sheet_network[strand_start]['bps']:

        # Its an parallel ladder if j + 1 is paired with strand_start + 1

        if self.sheet_network[j]['next'] != None \
           and  strand_start + 1 in self.sheet_network[j + 1]['bps']:
          num_parallel += 1

        else:
          num_antiparallel += 1

    if num_parallel == 0:
      self.type = 'antiparallel'
    elif num_antiparallel == 0:
      self.type = 'parallel'
    else:
      self.type = 'mixed'

  def parameterize(self):
    '''Parameterize the BetaSheet.
    The parameters are the 'bond' lengths, angles and torsions between Ca of each strand,
    as well as the vectors from each Ca to its beta pateners. The parameters are saved in
    the sheet_network.
    '''
    # Format the sheet_network

    for node in self.sheet_network:
      node['length'] = None
      node['angle'] = None
      node['torsion'] = None
      node['bp_vectors'] = []

    # Get the parameters
    
    i = 0
    while i < len(self.sheet_network):
      
      # Find the start and end of this strand

      strand_start = i
      strand_end = i

      while self.sheet_network[strand_end]['next']:
        strand_end += 1

      # Get bond lengths

      for j in range(strand_start, strand_end):
        self.sheet_network[j]['length'] = self.sheet_network[j + 1]['residue']['CA'] \
                - self.sheet_network[j]['residue']['CA']

      # Get bond angles

      for j in range(strand_start + 1, strand_end):
        self.sheet_network[j]['angle'] = np.degrees(PDB.calc_angle( 
                self.sheet_network[j - 1]['residue']['CA'].get_vector(),
                self.sheet_network[j]['residue']['CA'].get_vector(),
                self.sheet_network[j + 1]['residue']['CA'].get_vector()))

      # Get bond torsions

      for j in range(strand_start + 1, strand_end - 1):
        self.sheet_network[j]['torsion'] = np.degrees(PDB.calc_dihedral(
            self.sheet_network[j - 1]['residue']['CA'].get_vector(),
            self.sheet_network[j]['residue']['CA'].get_vector(),
            self.sheet_network[j + 1]['residue']['CA'].get_vector(),
            self.sheet_network[j + 2]['residue']['CA'].get_vector()))

      # Get the bp vectors

      for j in range(strand_start + 1, strand_end):
        stub = geometry.get_stub_matrix(self.sheet_network[j + 1]['residue']['CA'].get_coord(),
                self.sheet_network[j]['residue']['CA'].get_coord(),
                self.sheet_network[j - 1]['residue']['CA'].get_coord())

        for k in self.sheet_network[j]['bps']:
          bp_v_global = self.sheet_network[k]['residue']['CA'].get_coord() \
                  - self.sheet_network[j]['residue']['CA'].get_coord()

          self.sheet_network[j]['bp_vectors'].append(list(np.dot(np.array(stub.T), bp_v_global)))

      i = strand_end + 1


class Loop(SecondaryStructure):
  '''Class that represents a loop.'''
  def __init__(self, residue_list, ss_type):
    self.residue_list = residue_list
    self.type = ss_type


