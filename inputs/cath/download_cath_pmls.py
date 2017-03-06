#!/usr/bin/env python3

'''Download all superfamily structure superpositions of
all homologous superfamilies in the CATH database.
'''

import subprocess


if __name__ == '__main__':

  # Read the superfamily ids from a text file

  superfamily_ids = []

  with open('cath-superfamily-list.txt', 'r') as f:
    f.readline() # Ignore the first line

    while True:
      line = f.readline()
      if not line: break

      superfamily_ids.append(line.split()[0])

  # Download the pml files

  for sf_id in superfamily_ids:
    subprocess.call(['wget',
        'http://www.cathdb.info/version/v4_1_0/superfamily/{0}/superposition/pymol'.format(sf_id),
        '-O', 'cath-superfamily-superposition-{0}.pml'.format(sf_id)])
