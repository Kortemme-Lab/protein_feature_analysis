#!/usr/bin/env python3

import pytest

import numpy as np

from ProteinFeatureAnalyzer.features.geometry import *


def test_rotation_representation_conversion():

  rot_m = np.array([[ 0.25581,  -0.77351,   0.57986],
                    [-0.85333,  -0.46255,  -0.24057],
                    [ 0.45429,  -0.43327,  -0.77839]]) 

  # Test rotation matrix to euler angles

  theta_x, theta_y, theta_z = rotation_matrix_to_euler_angles(rot_m)
  
  assert(abs(-2.6337 - theta_x) < 0.001)
  assert(abs(-0.47158 - theta_y) < 0.001)
  assert(abs(-1.2795 - theta_z) < 0.001)

  # Test euler angles to rotation matrix

  rot_m2 = euler_angles_to_rotation_matrix(theta_x, theta_y, theta_z)

  for i in range(3):
    for j in range(3):
      assert(abs(rot_m[i][j] - rot_m2[i][j]) < 0.001)
