import os
import io
from datetime import timedelta

import numpy as np
import pandas as pd
import Bio.PDB as PDB
from flufl.lock import Lock


class Feature:
  '''Base class for features.'''

  def __init__(self):
    self.feature_list = []

  def list_my_jobs(self, input_path, total_num_threads, my_id):
    '''List all the inputs that should be hanled by the running thread.'''
    all_jobs = os.listdir(input_path)
    return [all_jobs[i] for i in range(len(all_jobs)) if i % total_num_threads == my_id]

  def append_to_csv(self, dataframe, file_name):
    '''Append a pandas dataframe to a csv file. This function is thread save.'''
    
    # Make a lock

    lock_name = os.path.join(file_name + '.lock')
    lock = Lock(lock_name)
    lock.lifetime = timedelta(minutes=10)

    # Write to the result file with a lock
    
    with lock:
      with open(file_name, 'a+') as f:
        dataframe.to_csv(f, header=False, index=False)
