import subprocess

import numpy as np
import Bio.PDB as PDB
import networkx as nx
import cylinder_fitting

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
    self.init_direction_vectors()
    self.init_angles()
    self.init_dihedrals()

  def init_direction_vectors(self):
    '''Initialize direction vectors corresponding to phi/psi torsions.'''
    
    # Get all peptide bond frames

    pp_bond_frames = []

    for i in range(len(self.residue_list) - 1):

        pp_bond_frames.append(np.array(geometry.get_stub_matrix(self.residue_list[i]['CA'].get_coord(),
            self.residue_list[i]['N'].get_coord(), self.residue_list[i + 1]['C'].get_coord())))

    # Get all transformations of peptide bonds

    pp_transformations = []

    for i in range(len(pp_bond_frames) - 1):

        pp_transformations.append(np.dot(pp_bond_frames[i + 1], np.transpose(pp_bond_frames[i])))

    # Get the direction vectors

    self.direction_vectors = [geometry.rotation_matrix_to_axis_and_angle(M)[0]
                                for M in pp_transformations]

  def init_angles(self):
    '''Initialize angles between direction vectors.'''
    self.angles = []

    for i in range(len(self.direction_vectors) - 1):
        self.angles.append(geometry.angle(self.direction_vectors[i], self.direction_vectors[i + 1]))

  def init_dihedrals(self):
    '''Initialize dihedrals among direction vectors'''
    self.dihedrals = []

    for i in range(len(self.direction_vectors) - 2):
        self.dihedrals.append(geometry.dihedral(np.zeros(3), np.zeros(3) + self.direction_vectors[i],
            np.zeros(3) + self.direction_vectors[i] + self.direction_vectors[i + 1],
            np.zeros(3) + self.direction_vectors[i] + self.direction_vectors[i + 1] + self.direction_vectors[i + 2]))

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
    self.init_gaussian_curvatures()
    self.fit_to_cylinder()

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
        
      self.graph.add_node(res, side=1000, terminal=1000, mismatch=1000, 
              gauss_cur='null', edge_cur1='null', edge_cur2='null', edge_cur3='null', edge_cur4='null',
              cylinder_radius='null')

    for strand in self.strand_list:
      for i in range(len(strand.residue_list) - 1):
        self.graph.add_edge(strand.residue_list[i],
                strand.residue_list[i + 1], edge_type='pp_bond')

    # Add the beta pairs

    for res in residues:
      d = dssp_dict[(res.get_parent().get_id(), res.get_id()[1])]
      
      for value in ['bp1', 'bp2']:
        if d[value] > 0:
          key = key_map[d[value]]
          res2 = model[key[0]][key[1]]

          if res2 in residues:
            self.graph.add_edge(res, res2, edge_type='bp')

    # Add the hydrogen bonds

    for res in residues:
      d = dssp_dict[(res.get_parent().get_id(), res.get_id()[1])]

      # Add a hydrogen bond if the H bond energy is below a limit

      HB_ENERGY_LIMIT = -1.5

      for value in ['nh_to_o1', 'nh_to_o2', 'o_to_nh1', 'o_to_nh2']:
        if d[value][1] < HB_ENERGY_LIMIT:
          key = key_map[d['seq_num'] + d[value][0]]
          res2 = model[key[0]][key[1]]
          
          if res2 in residues:
            self.graph.add_edge(res, res2, edge_type=value)
        
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
                    d['edge_type'] in ['nh_to_o1', 'nh_to_o2', 'o_to_nh1', 'o_to_nh2']]
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

  def init_gaussian_curvatures(self):
    '''Calculate gaussian curvatures of residues that
    are in the internal of a sheet.
    '''
    for res in self.graph.nodes():
      
      # Get the +2 residue

      p2_res = self.get_next_node(res)
      p2_res = self.get_next_node(p2_res) if p2_res else None
      if p2_res is None: continue

      # Get the -2 residue
      
      m2_res = self.get_prev_node(res)
      m2_res = self.get_prev_node(m2_res) if m2_res else None
      if m2_res is None: continue

      # Get the beta pair residues

      bp_res = set()

      for edge in self.graph.out_edges(res):
        if 'bp' in [d['edge_type'] for d in self.graph.get_edge_data(*edge).values()]:
          bp_res.add(edge[1])
          
      if len(bp_res) < 2:
        continue
        
      gaussian_curvature, edge_curvatures= geometry.calc_discrete_curvatures(
              res['CA'].get_coord(), [p2_res['CA'].get_coord(), bp_res.pop()['CA'].get_coord(),
                  m2_res['CA'].get_coord(), bp_res.pop()['CA'].get_coord()])

      self.graph.node[res]['gauss_cur'] = np.degrees(gaussian_curvature)
      self.graph.node[res]['edge_cur1'] = np.degrees(edge_curvatures[0])
      self.graph.node[res]['edge_cur2'] = np.degrees(edge_curvatures[1])
      self.graph.node[res]['edge_cur3'] = np.degrees(edge_curvatures[2])
      self.graph.node[res]['edge_cur4'] = np.degrees(edge_curvatures[3])

  def fit_to_cylinder(self):
    '''Fit local patches of beta sheets to a cylinder. Calculate
    the angle between the cylinder axies and the beta strand direction.
    Also calculate the local bending angle.
    '''
    for res in self.graph.nodes():

      # Get the +2 residue

      p2_res = self.get_next_node(res)
      p2_res = self.get_next_node(p2_res) if p2_res else None
      if p2_res is None: continue

      # Get the -2 residue
      
      m2_res = self.get_prev_node(res)
      m2_res = self.get_prev_node(m2_res) if m2_res else None
      if m2_res is None: continue

      # Get the beta pair residues for all 3 residues on the central strand

      bp_res = set()

      for central_res in [res, p2_res, m2_res]:
        for edge in self.graph.out_edges(central_res):
          if 'bp' in [d['edge_type'] for d in self.graph.get_edge_data(*edge).values()]:
            bp_res.add(edge[1])
     
      # At least 9 residues are neeeded for cylinder fitting

      if len(bp_res) < 6:
        continue

      # Get the beta pair residues for the centeral residue

      central_bp_res = set()

      for edge in self.graph.out_edges(res):
        if 'bp' in [d['edge_type'] for d in self.graph.get_edge_data(*edge).values()]:
          central_bp_res.add(edge[1])
     
      if len(central_bp_res) < 2:
        continue

      central_bp_res = list(central_bp_res)

      # Fit the CA atoms to a cylinder

      cas = [r['CA'].get_coord() for r in list(bp_res) + [res, p2_res, m2_res]]

      cas += [(res['CA'].get_coord() + p2_res['CA'].get_coord()) / 2, 
              (res['CA'].get_coord() + m2_res['CA'].get_coord()) / 2]

      w_fit, C_fit, r_fit, fit_err = cylinder_fitting.fit(cas)

      # Get the direction of the central strand

      strand_direction = geometry.normalize(p2_res['CA'].get_coord() - m2_res['CA'].get_coord())
      if np.dot(w_fit, strand_direction) < 0:
          strand_direction = -strand_direction

      # Get the sign of the angle between the central strand and the cylinder axis
      # positive for left handed folding and negative for right handed

      cr_v = geometry.projection(w_fit, res['CA'].get_coord() - C_fit)
      saa_sign = np.sign(np.dot(cr_v, np.cross(strand_direction, w_fit)))

      # Caculate the folding angle from the cylinder radius

      arc_strand = np.linalg.norm(geometry.projection(w_fit, 
                    p2_res['CA'].get_coord() - m2_res['CA'].get_coord()))
      arc_pair = np.linalg.norm(geometry.projection(w_fit, 
                    central_bp_res[0]['CA'].get_coord() -central_bp_res[1]['CA'].get_coord()))

      s_half = max(arc_strand, arc_pair) / (2 * r_fit)

      if s_half > 1: continue

      folding_angle = np.arcsin(s_half)
      
      #print("Fitting RMSD =", cylinder_fitting.fitting_rmsd(w_fit, C_fit, r_fit, cas)) ###DEBUG
      #print("Radius =", r_fit)###DEBUG
      #print("Arc =", max(arc_strand, arc_pair))
      #print("Folding_anlge =", folding_angle)###DEBUG
      #cylinder_fitting.show_fit(w_fit, C_fit, r_fit, cas) ###DEBUG

      self.graph.node[res]['cylinder_radius'] = r_fit
      self.graph.node[res]['cylinder_strand_angle'] = saa_sign * geometry.angle(w_fit, strand_direction)
      self.graph.node[res]['cylinder_fitting_rmsd'] = cylinder_fitting.fitting_rmsd(w_fit, C_fit, r_fit, cas)
      self.graph.node[res]['folding_angle'] = folding_angle

class Loop(SecondaryStructure):
  '''Class that represents a loop.'''
  def __init__(self, residue_list, ss_type):
    self.residue_list = residue_list
    self.type = ss_type

