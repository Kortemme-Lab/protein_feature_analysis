import os

import Bio.PDB as PDB

from .Feature import Feature


class StructuralHomologFeature(Feature):
  '''Analyze features of structrual homologs.'''
  
  def __init__(self):
    super().__init__()
    self.proteins = []

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
          os.mkdir(scratch_path, exist_ok=True)

        # Extract features from one file  

        self.extract_from_one_file(os.path.join(input_path, f), scratch_path)

  def extract_from_one_file(self, pml_file, scratch_path):
    '''Extract structrual homolog features from a .pml file of superposed
    homologous superfamilies from the CATH database.
    '''
    self.proteins = []
    
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

            pdb_lines.append(line.strip().strip('\\'))
         
          # Make a pdb file of the structure for DSSP analysis
          
          with open(os.path.join(scratch_path, pdb_id + '.pdb'), 'w') as pdb_f:
            pdb_f.write('\n'.join(pdb_lines))

          self.proteins.append(
                {'structure':self.structure_from_pdb_string('\n'.join(pdb_lines), pdb_id)} )

    for p in self.proteins:
      self.find_secondary_structures(p)

  def find_secondary_structures(self, protein_dict):
    '''Find secondary structures of a protein.
    
    Arguements:
     - protein_dict - a dictionary to store informations of a protein
    '''
    pass
