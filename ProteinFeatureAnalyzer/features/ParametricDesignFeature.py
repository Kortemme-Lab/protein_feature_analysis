import os
import json

import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import pandas as pd
import Bio.PDB as PDB

from .Feature import Feature
from . import topology
from . import geometry
from . import secondary_structures
from . import data_loading


class ParametricDesignFeature(Feature):
  '''Analyze features for parametric design.'''
  
  def __init__(self, dssp_path):
    super().__init__()
    self.superfamilies = []
    self.dssp_path = dssp_path

  def extract(self, input_path, output_path, total_num_threads=1, my_id=0):
    '''Extract structrual homolog features from structures in the input path.
    The input data should be stored in .pml files of superposed homologous 
    superfamilies from the CATH database.
    '''

    self.superfamilies = data_loading.load_data_from_cath_pmls(input_path, output_path,
                            self.list_my_jobs(input_path, total_num_threads, my_id), self.dssp_path)


