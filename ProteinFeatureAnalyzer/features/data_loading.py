import os
import io

import Bio.PDB as PDB

from . import topology
from . import secondary_structures


def structure_from_pdb_file(file_path, name=''):
  '''Read the structure stored in a PDB file.'''
  parser = PDB.PDBParser()
  return parser.get_structure(name, file_path)

def structure_from_pdb_string(pdb_string, name=''):
  '''Read the structure stored in a PDB string.'''
  parser = PDB.PDBParser()
  pdb_sf = io.StringIO(pdb_string)
  return parser.get_structure(name, pdb_sf)


def load_data_from_cath_pmls(input_path, output_path, job_list, dssp_path):
  '''Load data from structures in the input path.
  The input data should be stored in .pml files of superposed homologous 
  superfamilies from the CATH database.
  '''

  superfamilies = []

  for f in job_list:
    if f.endswith('.pml'):
    
      # Make a scratch directory

      scratch_path = os.path.join(output_path, f[0:-4])
      if not os.path.exists(scratch_path):
        os.mkdir(scratch_path)

      # Load data from one file  

      load_from_one_cath_pml_file(os.path.join(input_path, f), scratch_path, superfamilies, dssp_path)

  return superfamilies

def load_from_one_cath_pml_file(pml_file, scratch_path, superfamilies, dssp_path):
  '''Load data from a .pml file of superposed
  homologous superfamilies from the CATH database.
  '''
  superfamilies.append([])
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
       
        structure = structure_from_pdb_string('\n'.join(pdb_lines), pdb_id)

        # Store structures without chain breaks

        if len(topology.find_structure_chain_breaks(structure)) == 0:
          structure_path = os.path.join(scratch_path, pdb_id + '.pdb')

          io = PDB.PDBIO()
          io.set_structure(structure)
          io.save(structure_path)

          candidate_proteins.append({'structure' : structure, 'path' : structure_path})

  for p in candidate_proteins:
    try:  
      find_secondary_structures(p, dssp_path)
    except:
      continue
    superfamilies[-1].append(p) # Add a protein to a superfamily if there's no exception

def find_secondary_structures(protein_dict, dssp_path):
  '''Find secondary structures of a protein.
  
  Arguements:
   - protein_dict - a dictionary to store informations of a protein
  '''
  protein_dict['dssp_dict'], protein_dict['dssp_key_map'] = \
    secondary_structures.make_dssp_dict(protein_dict['path'], dssp_path)

  protein_dict['ss_list'], protein_dict['sheet_list'] = \
    secondary_structures.pack_dssp_dict_into_ss_list(protein_dict['structure'][0],
          protein_dict['dssp_dict'], protein_dict['dssp_key_map'])
