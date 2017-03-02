import os

import Bio.PDB as PDB

from .Feature import Feature


class StructuralHomologFeature(Feature):
  '''Analyze features of structrual homologs.'''
  
  def __init__(self):
    super().__init__()
    self.structures = []

  def extract(self, input_path, total_num_threads=1, my_id=0):
    '''Extract structrual homolog features from structures in the input path.
    The input data should be stored in .pml files of superposed homologous 
    superfamilies from the CATH database.
    '''
    for f in self.list_my_jobs(input_path, total_num_threads, my_id):
      if f.endswith('.pml'):
        self.extract_from_one_file(os.path.join(input_path, f))

  def extract_from_one_file(self, pml_file):
    '''Extract structrual homolog features from a .pml file of superposed
    homologous superfamilies from the CATH database.
    '''
    self.structures = []
    
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
              pdb_id = line.strip()[5:13]
              break

            pdb_lines.append(line.strip().strip('\\'))
          
          self.structures.append( self.structure_from_pdb_string('\n'.join(pdb_lines), pdb_id) )

