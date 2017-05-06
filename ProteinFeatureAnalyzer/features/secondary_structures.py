import subprocess

import numpy as np
import Bio.PDB as PDB
import networkx as nx

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
    nh_to_o1 = l[39:50].split(',')     # Hydrogen bond from NH to O 1
    o_to_nh1 = l[50:61].split(',')     # Hydrogen bond from O to NH 1
    nh_to_o2 = l[61:72].split(',')     # Hydrogen bond from NH to O 2
    o_to_nh2 = l[72:83].split(',')     # Hydrogen bond from O to NH 2

    key_map[seq_num] = (chain_id, res_id)
    dssp_dict[(chain_id, res_id)] = {'seq_num':seq_num, 'ss':ss,
            'bbl1':bbl1, 'bbl2':bbl2, 'bp1':bp1, 'bp2':bp2, 'bsl':bsl,
            'nh_to_o1':(int(nh_to_o1[0]), float(nh_to_o1[1])),
            'o_to_nh1':(int(o_to_nh1[0]), float(o_to_nh1[1])),
            'nh_to_o2':(int(nh_to_o2[0]), float(nh_to_o2[1])),
            'o_to_nh2':(int(o_to_nh2[0]), float(o_to_nh2[1])) }

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
        dssp_dict, key_map, model))

  return ss_list, sheet_list

class SecondaryStructure:
  '''Base class for secondary structures.'''
  def __init__(self):
    pass

class AlphaHelix(SecondaryStructure):
  '''Class that represents an alpha helix.'''
  def __init__(self, residue_list):
    self.residue_list = residue_list

class BetaStrand(SecondaryStructure):
  '''Class that represents a beta strand.'''
  def __init__(self, residue_list, sheet_id):
    self.residue_list = residue_list
    self.sheet_id = sheet_id

class BetaSheet(SecondaryStructure):
  '''Class that represents a beta sheet.'''
  def __init__(self, sheet_id, strand_list, dssp_dict, key_map, model):
    self.sheet_id = sheet_id
    self.strand_list = strand_list
    self.init_sheet_graph(dssp_dict, key_map, model)
    self.init_boundaries(dssp_dict)
    self.init_mismatches()
    self.propagate_node_attributes()
    self.determine_sheet_type()

    for res in self.graph.nodes(): ###DEBUG
      print(self.graph.node[res])

  def init_sheet_graph(self, dssp_dict, key_map, model):
    '''Initialize a graph to represent a beta sheet.
    Nodes are residues and edges are peptide bonds or
    hydrogen bonds.
    '''
    self.graph = nx.MultiDiGraph()

    # Add the strands
    
    residues = [res for strand in self.strand_list for res in strand.residue_list]
    for res in residues:
      
      # The attributes side, terminal and mismatch tell the shortest distance
      # of a residue to side, terminal and mismatch residues.
        
      self.graph.add_node(res, side=1000, terminal=1000, mismatch=1000)

    for strand in self.strand_list:
      for i in range(len(strand.residue_list) - 1):
        self.graph.add_edge(strand.residue_list[i],
                strand.residue_list[i + 1], edge_type='pp_bond')

    # Add the hydrogen bonds

    for res in residues:
      d = dssp_dict[(res.get_parent().get_id(), res.get_id()[1])]

      # Add a hydrogen bond if the H bond energy is below a limit

      HB_ENERGY_LIMIT = -1.5

      if d['nh_to_o1'][1] < HB_ENERGY_LIMIT:
        key = key_map[d['seq_num'] + d['nh_to_o1'][0]]
        res2 = model[key[0]][key[1]]
        
        if res2 in residues:
          self.graph.add_edge(res, res2, edge_type='nh_to_o1')
        
      if d['o_to_nh1'][1] < HB_ENERGY_LIMIT:
        key = key_map[d['seq_num'] + d['o_to_nh1'][0]]
        res2 = model[key[0]][key[1]]
        
        if res2 in residues:
          self.graph.add_edge(res, res2, edge_type='o_to_nh1')
        
      if d['nh_to_o2'][1] < HB_ENERGY_LIMIT:
        key = key_map[d['seq_num'] + d['nh_to_o2'][0]]
        res2 = model[key[0]][key[1]]
        
        if res2 in residues:
          self.graph.add_edge(res, res2, edge_type='nh_to_o2')
        
      if d['o_to_nh2'][1] < HB_ENERGY_LIMIT:
        key = key_map[d['seq_num'] + d['o_to_nh2'][0]]
        res2 = model[key[0]][key[1]]
        
        if res2 in residues:
          self.graph.add_edge(res, res2, edge_type='o_to_nh2')

  def get_prev_node(self, residue):
    '''Get the previous node of a residue in the beta sheet graph.
    Return None in the previous node doesn't exist.
    '''
    for edge in self.graph.in_edges(residue):
      edge_types = [d['edge_type'] for d in self.graph.get_edge_data(*edge).values()]

      if 'pp_bond' in edge_types:
        return edge[0]

    return None

  def get_next_node(self, residue):
    '''Get the next node of a residue in the beta sheet graph.
    Return None if next node doesn't exist.
    '''
    for edge in self.graph.out_edges(residue): 
      edge_types = [d['edge_type'] for d in self.graph.get_edge_data(*edge).values()]
      
      if 'pp_bond' in edge_types:
        return edge[1]

    return None

  def init_boundaries(self, dssp_dict): 
    '''Initialize the boundary residues of the graph.'''

    for res in self.graph.nodes():
      key = (res.get_parent().get_id(), res.get_id()[1])
      
      # Set side residues which has only one beta sheet pair residue

      if dssp_dict[key]['bp1'] == 0 or dssp_dict[key]['bp2'] == 0:
        self.graph.node[res]['side'] = 0

      # Set terminal residues

      if self.get_next_node(res) is None or self.get_prev_node(res) is None:
        self.graph.node[res]['terminal'] = 0

  def num_hbonds(self, residue):
    '''Return the number of Hbonds of a residue.'''
    hbonds = []
    
    for edge in self.graph.out_edges(residue):
      hbonds += [d for d in self.graph.get_edge_data(*edge).values() if
                    d['edge_type'] != 'pp_bond']
    return len(hbonds)
  
  def init_mismatches(self):
    '''Initialize the mismatched residues.
    Internal and terminal residues without both N and O Hbonded are regarded as mismatches.
    Side residue are mismatches unless they have at leat one Hbond or both its previous
    and next residues have two Hbonds.
    '''
    for res in self.graph.nodes():

      if self.graph.node[res]['side'] > 0 and self.num_hbonds(res) < 2:
        self.graph.node[res]['mismatch'] = 0

      elif self.graph.node[res]['side'] == 0 and self.num_hbonds(res) < 1:
        prev_r = self.get_prev_node(res)
        next_r = self.get_next_node(res)

        if ((prev_r is not None) and self.num_hbonds(prev_r) < 2) \
                or ((next_r is not None) and self.num_hbonds(next_r) < 2):
          self.graph.node[res]['mismatch'] = 0

  def propagate_node_attributes(self):
    '''Propagate the the side, teminal and mismatch attributes throughout
    the whole graph.
    '''
    def propagete_recursively(node):
      '''Propagate the attributes recursively.'''
      updated_nodes = set()

      for neighbor in self.graph.neighbors(node):
        
        for attribute in ['side', 'terminal', 'mismatch']:
          if self.graph.node[neighbor][attribute] > self.graph.node[node][attribute] + 1:
            self.graph.node[neighbor][attribute] = self.graph.node[node][attribute] + 1
            updated_nodes.add(neighbor)

      for neighbor in updated_nodes:
        propagete_recursively(neighbor)

    for node in self.graph.nodes():
      propagete_recursively(node)

  def determine_sheet_type(self):
    '''Determine the type of the beta sheet.'''
    directions = [] # Save all the directions found

    for node1 in self.graph.nodes():
      next1 = self.get_next_node(node1)

      for edge in self.graph.out_edges(node1):
        edge_types = [d['edge_type'] for d in self.graph.get_edge_data(*edge).values()]
        
        if 'pp_bond' in edge_types:
          continue
        
        node2 = edge[1]
        next2 = self.get_next_node(node2)

        if bool(next1 and next2):
          v1 = next1['CA'].get_coord() - node1['CA'].get_coord()
          v2 = next2['CA'].get_coord() - node2['CA'].get_coord()

          directions.append(np.dot(v1, v2) > 0)

    if not (True in directions):
      self.type = 'antiparallel'
    elif not (False in directions):
      self.type = 'parallel'
    else:
      self.type = 'mixed'

class Loop(SecondaryStructure):
  '''Class that represents a loop.'''
  def __init__(self, residue_list, ss_type):
    self.residue_list = residue_list
    self.type = ss_type

