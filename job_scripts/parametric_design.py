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

  feature = PFA.features.ParametricDesignFeature(dssp_path)

  feature.extract(input_path, data_path, num_jobs, job_id)
  
  #feature.save_alpha_helix_features(data_path)
  feature.save_beta_sheet_features(data_path)
  
  #feature.load_alpha_helix_features(data_path)
  #feature.visualize_alpha_helix_features('angles')
  #feature.load_beta_sheet_features(data_path)
  #feature.visualize_beta_sheet_features('cylinder_radius')
  
