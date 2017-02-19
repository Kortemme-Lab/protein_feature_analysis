'''Helper functions for machine learning.'''

import numpy as np


def angle_to_cos_sin(angle):
  '''Convert an angle in radian to its cos and sin.
  This useful to preserve the topology of periodic features.
  '''
  return [np.cos(angle), np.sin(angle)]

def cos_sin_to_angle(cos, sin):
  '''Get an angle from its cos and sin.
  The range of the angle is [-pi, -pi]
  '''
  return np.arctan2(sin, cos)
