#!/usr/bin/env python3

import sys

import ProteinFeatureAnalyzer as PFA

if __name__ == '__main__':
  data_path = sys.argv[1]
  input_path = sys.argv[2]

  feature = PFA.features.RamachandranFeature()
  feature.extract(input_path)
  feature.visualize()
