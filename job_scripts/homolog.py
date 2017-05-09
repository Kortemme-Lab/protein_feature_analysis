#!/usr/bin/env python3

import sys

import ProteinFeatureAnalyzer as PFA

if __name__ == '__main__':
  
  data_path = sys.argv[1]
  input_path = sys.argv[2]
  
  num_jobs = 1
  job_id = 0
  if len(sys.argv) > 4:
    num_jobs = int(sys.argv[3])
    job_id = int(sys.argv[4]) - 1
  
  dssp_path = './dependencies/dssp-2.0.4-linux-i386'

  feature = PFA.features.StructuralHomologFeature(dssp_path)
  
  feature.extract(input_path, data_path, num_jobs, job_id)
  
  #feature.save_ss_composition_features(data_path)
  #feature.save_ss_sizes_features(data_path)
  #feature.save_beta_sheet_features(data_path)

  #feature.load_ss_composition_features(data_path)
  #feature.visualize_ss_composition_features('num_alpha_helix')
  #feature.load_ss_sizes_features(data_path)
  #feature.visualize_ss_sizes_features(ss='beta_sheet', s_type='antiparallel', fig_save_path=data_path)
  feature.load_beta_sheet_features(data_path)

  #for st in ['parallel', 'antiparallel']:
  #  for dof in ['phi', 'psi']:
  #    for t in [[(0, 10000), (0, 10000), (0, 10000)], [(0, 10000), (1, 10000), (1, 10000)], [(1, 10000), (0, 10000), (1, 10000)], [(1, 10000), (1, 10000), (0, 10000)], [(1, 10000), (1, 10000), (1, 10000)]]:
  #      print('<td>{0}</td><td>{1}</td><td></td><td></td><td></td><td>{2}</td><td>{3}</td><td>{4}</td>'.format(st, dof, t[0], t[1], t[2]))
  #      feature.visualize_beta_sheet_features(st, dof, side_threshold=t[0], terminal_threshold=t[1], mismatch_threshold=t[2], fig_save_path=data_path)
  
  feature.visualize_beta_sheet_features('parallel', 'edge_cur', side_threshold=(0,10000), terminal_threshold=(0, 10000), mismatch_threshold=(0, 10000))
  #feature.visualize_beta_sheet_features('antiparallel', 'phi', feature2='psi', side_threshold=(0,10000), terminal_threshold=(0, 10000), mismatch_threshold=(0, 10000))
  
