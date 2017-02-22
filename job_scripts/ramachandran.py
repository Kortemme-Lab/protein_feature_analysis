#!/usr/bin/env python3

import sys

import ProteinFeatureAnalyzer as PFA

if __name__ == '__main__':
  data_path = sys.argv[1]
  input_path = sys.argv[2]

  feature = PFA.features.RamachandranFeature()
  #feature.extract(input_path)
  #feature.save(data_path)
  feature.load(data_path)
  feature.learn(transform_features=True)
  #feature.learn("IsolationForest", transform_features=True)
  feature.calculate_space_reduction(transform_features=True)
  feature.visualize(transform_features=True)
