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
  
  #feature.extract(input_path, data_path, num_jobs, job_id)
  
  #feature.save_ss_composition_features(data_path)
  #feature.save_ss_sizes_features(data_path)
  #feature.save_alpha_helix_parameterization_features(data_path)

  #feature.load_ss_composition_features(data_path)
  #feature.visualize_ss_composition_features('num_alpha_helix')
  #feature.load_ss_sizes_features(data_path)
  #feature.visualize_ss_sizes_features(ss='loop', s_type='beta-beta', fig_save_path=data_path)
  feature.load_alpha_helix_parameterization_features(data_path)
  feature.visualize_alpha_helix_parameterization_features(v_type='a_t', position_shift=5)
  
