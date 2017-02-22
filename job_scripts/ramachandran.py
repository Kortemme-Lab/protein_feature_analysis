#!/usr/bin/env python3

import sys

import ProteinFeatureAnalyzer as PFA

if __name__ == '__main__':
  data_path = sys.argv[1]
  input_path = sys.argv[2]

  tf = True

  feature = PFA.features.RamachandranFeature()
  #feature.extract(input_path)
  #feature.save(data_path)
  feature.load(data_path)
  
  #feature.learn(transform_features=tf)
  #feature.learn("IsolationForest", transform_features=tf)
  #feature.calculate_space_reduction(transform_features=tf)
  
  feature.density_estimate()
  
  feature.visualize(transform_features=tf)
