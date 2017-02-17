import numpy as np
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
