#!/usr/bin/env python3

import sys

import ProteinFeatureAnalyzer as PFA

if __name__ == '__main__':
  data_path = sys.argv[1]
  input_path = sys.argv[2]
  
  dssp_path = './dependencies/dssp-2.0.4-linux-i386'

  feature = PFA.features.StructuralHomologFeature(dssp_path)
  
  #feature.extract(input_path, data_path)

  #feature.save_ss_composition_features(data_path)
  feature.load_ss_composition_features(data_path)
  
