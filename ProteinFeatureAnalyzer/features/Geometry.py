import numpy as np
import scipy
import scipy.spatial
import Bio.PDB as PDB


def get_phi(chain, residue):
  '''Calculate the phi torsion of a residue.'''
  
  # Get the previous residue

  res_id = residue.get_id()
  prev_res = chain[res_id[1] - 1]
  
  prev_flag = prev_res.get_id()[0]
  if prev_flag == 'W' or prev_flag.startswith('H_'):
    raise Exception('Hetero residue type!')
 
  # Calculate the torsion

  c_prev = prev_res['C'].get_vector()
  n = residue['N'].get_vector()
  ca = residue['CA'].get_vector()
  c = residue['C'].get_vector()

  return PDB.calc_dihedral(c_prev, n, ca, c) 

def get_psi(chain, residue):
  '''Calculate the psi torsion of a residue.'''
  
  # Get the next residue

  res_id = residue.get_id()
  next_res = chain[res_id[1] + 1]
  
  next_flag = next_res.get_id()[0]
  if next_flag == 'W' or next_flag.startswith('H_'):
    raise Exception('Hetero residue type!')
  
  # Calculate the torsion

  n = residue['N'].get_vector()
  ca = residue['CA'].get_vector()
  c = residue['C'].get_vector()
  n_next = next_res['N'].get_vector()

  return PDB.calc_dihedral(n, ca, c, n_next)

def get_distance_matrix(atom_list):
  '''Get the distance matrix of a list of atoms.'''
  return scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(
      np.array([a.get_coord() for a in atom_list]), 'euclidean'))

def get_nearest_nonbonded_residues(structure):
  '''Return a list of 2-tuples. Each the first element of each tuple is a
  residue and the second element is its nearest nonbonded residue.'''
 
  # Get all CA atoms which are used to be the center of residues

  ca_list = []
  
  for residue in structure.get_residues():
    flag = residue.get_id()[0]
    if flag == 'W' or flag.startswith('H_'): continue
    
    for a in residue:
      if a.get_id() == 'CA':
        ca_list.append(a)

  ca_coords = [a.get_coord() for a in ca_list]

  # Make a KDTree for neighbor searching

  kd_tree = scipy.spatial.KDTree(ca_coords)

  # Find the nearest nonbonded neighbor of all residues

  nearest_nb_list = []

  for i in range(len(ca_list)):
    res1 = ca_list[i].get_parent()
    distance, indices = kd_tree.query(ca_coords[i], k=4)

    for j in range(1, 4):
      res2 = ca_list[indices[j]].get_parent() 

      if res1.get_id()[1] + 1 == res2.get_id()[1] \
         or res1.get_id()[1] - 1 == res2.get_id()[1] : # Bonded residues
           continue
    
    nearest_nb_list.append((res1, res2))   

  return nearest_nb_list


