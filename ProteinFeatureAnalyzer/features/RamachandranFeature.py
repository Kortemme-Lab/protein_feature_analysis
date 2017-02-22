import os

import matplotlib
matplotlib.use('TkAgg') 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Bio.PDB as PDB
from sklearn import svm
from sklearn.ensemble import IsolationForest

from .Feature import Feature
from . import geometry
from . import machine_learning

class RamachandranFeature(Feature):
  '''Analyze the phi/psi torsion distributions of proteins.'''
  
  def __init__(self):
    super().__init__()
    self.clf = None

  def extract(self, input_path, total_num_threads=1, my_id=0):
    '''Extract phi, psi angles from structures in the input path.'''
    for f in self.list_my_jobs(input_path, total_num_threads, my_id):
      if f.endswith('.pdb'):
        self.extract_from_one_file(os.path.join(input_path, f))

  def extract_from_one_file(self, pdb_file):
    '''Extract phi, psi angles from a pdb_file.'''
    structure = self.structure_from_pdb_file(pdb_file)

    for model in structure:
      for chain in model:
        for residue in chain:
          try:
            feature_dict = {'phi' : geometry.get_phi(chain, residue),
                            'psi' : geometry.get_psi(chain, residue)}
            self.feature_list.append(feature_dict)
          except:
            pass

  def visualize(self):
    '''Visualize the feature statistics.'''
    phis = [ d['phi'] for d in self.feature_list ]
    psis = [ d['psi'] for d in self.feature_list ]

    # Draw the decision boundary from the machine learning classifier
    
    if self.clf:
      xx, yy = np.meshgrid(np.linspace(-np.pi, np.pi, 200), np.linspace(-np.pi, np.pi, 200))
      transformed_data = [ machine_learning.angle_to_cos_sin(d[0]) + machine_learning.angle_to_cos_sin(d[1])
		for d in np.c_[xx.ravel(), yy.ravel()] ]

      Z = self.clf.decision_function(transformed_data)
      Z = Z.reshape(xx.shape)

      Z_pred = self.clf.predict(transformed_data)
      Z_pred = Z_pred.reshape(xx.shape)
      
      #plt.contour(xx, yy, Z, levels=[0], linewidths=2, colors='darkred')
      plt.contourf(xx, yy, Z, levels=np.linspace(Z.min(), Z.max(), 7), cmap=plt.cm.Blues_r)
      #plt.contourf(xx, yy, Z, levels=np.linspace(0, Z.max()), colors='orange')
      plt.contourf(xx, yy, Z_pred, levels=[0.9, 1.1], colors='orange')

    # Draw the data

    plt.scatter(phis, psis, c='green', s=5)

    # Plot the support vectors if the classifier is SVM

    if isinstance(self.clf, svm.OneClassSVM):
      s_phis = [ machine_learning.cos_sin_to_angle(v[0], v[1]) for v in self.clf.support_vectors_ ]
      s_psis = [ machine_learning.cos_sin_to_angle(v[2], v[3]) for v in self.clf.support_vectors_ ]
      plt.scatter(s_phis, s_psis, c='red')

    plt.axis([- np.pi, np.pi, - np.pi, np.pi])
    plt.show()

  def save(self, data_path):
    '''Save the data into a csv file.'''
    data = [ (d['phi'], d['psi']) for d in self.feature_list ]
    df = pd.DataFrame(data=data, columns=['phi', 'psi'])
    
    self.append_to_csv(df, os.path.join(data_path, 'rama_features.csv'))

  def load(self, data_path):
    '''Load data from a csv file.'''
    df = pd.read_csv(os.path.join(data_path, 'rama_features.csv'), header=None)
    
    for index, row in df.iterrows():
      self.feature_list.append({'phi':row[0], 'psi':row[1]})

  def learn(self, clf_type="OneClassSVM"):
    '''Learn the distribution with a machine learning classifier'''
    # Prepare the training data
    
    all_data = [machine_learning.angle_to_cos_sin(d['phi']) + machine_learning.angle_to_cos_sin(d['psi'])
            for d in self.feature_list]
    n_data = len(all_data)

    training_data = all_data[0:int(0.6 * n_data)]
    test_data = all_data[int(0.6 * n_data):int(0.8 * n_data)]
    cv_data = all_data[int(0.8 * n_data):n_data]

    # Train the classifier 

    if clf_type == "OneClassSVM":
      nus = [0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
      least_error = len(test_data)

      for i in range(len(nus)):
        print("nu = {0}".format(nus[i]))

        clf = svm.OneClassSVM(nu=nus[i], kernel="rbf", gamma='auto')
        clf.fit(training_data)
        
        predictions = clf.predict(training_data)
        print("{0}/{1} training error.".format(len(predictions[-1 == predictions]), len(training_data)))
        
        predictions = clf.predict(test_data)
        print("{0}/{1} test error.\n".format(len(predictions[-1 == predictions]), len(test_data)))

        if len(predictions[-1 == predictions]) < least_error:
          least_error = len(predictions[-1 == predictions])
          self.clf = clf
    
    elif clf_type == "IsolationForest": 
      self.clf = IsolationForest(max_samples=200,
			contamination=0.05, random_state=np.random.RandomState(42))
      self.clf.fit(training_data)
   
    # Print Training results
    
    predictions = self.clf.predict(cv_data)
    print("{0}/{1} cross validation error.".format(len(predictions[-1 == predictions]), len(cv_data)))

    if clf_type == "OneClassSVM":
      print("{0} support vectors found.".format(len(self.clf.support_)))

  def predict(self, input_data):
    '''Make a prediction for the input data with the machine learning classifier.
    input_data is a list of phi, psi angles.
    '''
    transformed_data = [machine_learning.angle_to_cos_sin(d[0]) + machine_learning.angle_to_cos_sin(d[1])
            for d in input_data]
    return self.clf.predict(transformed_data)

