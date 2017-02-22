#!/usr/bin/env python3

import sys

import numpy as np

import ProteinFeatureAnalyzer as PFA

if __name__ == '__main__':
  data_path = sys.argv[1]
  input_path = sys.argv[2]

  feature = PFA.features.BackboneMicroEnvironmentFeature()
  #feature.extract(input_path)
  #feature.save(data_path)
  feature.load(data_path)

  ######### MACHINE LEARNING #########
 
  feature.learn()
  
  ######### VISUALIZATION #########

  #feature.visualize()
  #feature.plot_shifts()
  #feature.plot_shift_length_histogram()
  #feature.scatter_plot_two_features(lambda d : d['phi1'], lambda d : d['psi1'], [-np.pi, np.pi, -np.pi, np.pi])
  #feature.scatter_plot_two_features(lambda d : d['phi2'], lambda d : d['psi2'], [-np.pi, np.pi, -np.pi, np.pi])
  
  #feature.scatter_plot_two_features(lambda d : d['phi1'], lambda d : d['phi2'], [-np.pi, np.pi, -np.pi, np.pi])
  #feature.scatter_plot_two_features(lambda d : d['psi1'], lambda d : d['psi2'], [-np.pi, np.pi, -np.pi, np.pi])
  
  #feature.scatter_plot_two_features(lambda d : d['phi1'], lambda d : d['psi2'], [-np.pi, np.pi, -np.pi, np.pi])
  #feature.scatter_plot_two_features(lambda d : d['psi1'], lambda d : d['phi2'], [-np.pi, np.pi, -np.pi, np.pi])
  
  #feature.scatter_plot_two_features(lambda d : d['shift'][0], lambda d : d['shift'][1])
  #feature.scatter_plot_two_features(lambda d : d['shift'][0], lambda d : d['shift'][2])
  #feature.scatter_plot_two_features(lambda d : d['shift'][1], lambda d : d['shift'][2])
 
  #feature.scatter_plot_two_features(lambda d : d['shift'][0], lambda d : d['theta_x'])
  #feature.scatter_plot_two_features(lambda d : d['shift'][1], lambda d : d['theta_y'])
  #feature.scatter_plot_two_features(lambda d : d['shift'][2], lambda d : d['theta_z'])

  tx_r = [-np.pi, np.pi] # Range of theta_x
  ty_r = [-np.pi / 2, np.pi / 2] # Range of theta_y
  tz_r = [-np.pi, np.pi] # Range of theta_z

  #feature.scatter_plot_two_features(lambda d : d['theta_x'], lambda d : d['theta_y'], tx_r + ty_r)
  #feature.scatter_plot_two_features(lambda d : d['theta_x'], lambda d : d['theta_z'], tx_r + tz_r)
  #feature.scatter_plot_two_features(lambda d : d['theta_y'], lambda d : d['theta_z'], ty_r + tz_r)
  
  #feature.scatter_plot_three_features(lambda d : d['theta_x'], lambda d : d['theta_y'], lambda d : d['theta_z'], tx_r + ty_r + tz_r)
