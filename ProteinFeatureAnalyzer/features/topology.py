import Bio.PDB as PDB


def find_structure_chain_breaks(structure):
  '''Find the chain breaks of a structure.
  A chain break is represented by a tuple (id_start, id_end),
  where from id_start + 1 to id_end - 1 are missing residues.
  '''
  chain_breaks = []

  for model in structure:
    for chain in model:

      residue_list = [r for r in chain.get_residues()]

      for i in range(len(residue_list) - 1):
        if residue_list[i].get_id()[1] + 1 != residue_list[i + 1].get_id()[1]:
          chain_breaks.append(((model.get_id(), chain.get_id(), residue_list[i].get_id()),
                               (model.get_id(), chain.get_id(), residue_list[i + 1].get_id())))

  return chain_breaks
