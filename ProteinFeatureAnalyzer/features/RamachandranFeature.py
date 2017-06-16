import os

import matplotlib
matplotlib.use('TkAgg') 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Bio.PDB as PDB
from sklearn import svm
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from sklearn import mixture

from .Feature import Feature
from . import data_loading
from . import geometry
from . import machine_learning

class RamachandranFeature(Feature):
  '''Analyze the phi/psi torsion distributions of proteins.'''
  
  def __init__(self):
    super().__init__()
    self.clf = None
    self.de = None

  def extract(self, input_path, total_num_threads=1, my_id=0):
    '''Extract phi, psi angles from structures in the input path.'''
    for f in self.list_my_jobs(input_path, total_num_threads, my_id):
      if f.endswith('.pdb'):
        self.extract_from_one_file(os.path.join(input_path, f))

  def extract_from_one_file(self, pdb_file):
    '''Extract phi, psi angles from a pdb_file.'''
    structure = data_loading.structure_from_pdb_file(pdb_file)

    for model in structure:
      for chain in model:
        for residue in chain:
          try:
            feature_dict = {'phi' : geometry.get_phi(chain, residue),
                            'psi' : geometry.get_psi(chain, residue)}
            self.feature_list.append(feature_dict)
          except:
            pass

  def visualize(self, transform_features=True):
    '''Visualize the feature statistics.'''
    phis = [ d['phi'] for d in self.feature_list ]
    psis = [ d['psi'] for d in self.feature_list ]

    # Prepare grid points

    xx, yy = np.meshgrid(np.linspace(-np.pi, np.pi, 200), np.linspace(-np.pi, np.pi, 200))
    
    transformed_data = np.c_[xx.ravel(), yy.ravel()]

    if transform_features:
      transformed_data = self.transform_features(transformed_data)
    
    # Draw the decision boundary from the machine learning classifier
    
    if self.clf:

      Z = self.clf.decision_function(transformed_data)
      Z = Z.reshape(xx.shape)

      Z_pred = self.clf.predict(transformed_data)
      Z_pred = Z_pred.reshape(xx.shape)
      
      #plt.contour(xx, yy, Z, levels=[0], linewidths=2, colors='darkred')
      plt.contourf(xx, yy, Z, levels=np.linspace(Z.min(), Z.max(), 20), cmap=plt.cm.Blues_r)
      #plt.contourf(xx, yy, Z, levels=np.linspace(0, Z.max()), colors='orange')
      plt.contourf(xx, yy, Z_pred, levels=[0.9, 1.1], colors='orange')

    # Draw the density estimation

    if self.de:
      
      Z = self.de.score_samples(transformed_data)
      Z = Z.reshape(xx.shape)
      plt.contourf(xx, yy, Z, levels=np.linspace(Z.min(), Z.max(), 7), cmap=plt.cm.Blues_r)

    # Draw the data

    plt.scatter(phis, psis, c='green', s=5)

    # Plot the support vectors if the classifier is SVM

    if isinstance(self.clf, svm.OneClassSVM):

      if transform_features:  
        s_phis = [ machine_learning.cos_sin_to_angle(v[0], v[1]) for v in self.clf.support_vectors_ ]
        s_psis = [ machine_learning.cos_sin_to_angle(v[2], v[3]) for v in self.clf.support_vectors_ ]
        plt.scatter(s_phis, s_psis, c='red')
      else:
        plt.scatter(self.clf.support_vectors_[:][0], self.clf.support_vectors_[:][1], c='red')

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

  def transform_features(self, feature_list):
    '''Transform feature representations. The arguement feature_list
    could be a list of dictionary or a list of list.
    '''
    if isinstance(feature_list[0], dict):
      return [machine_learning.angle_to_cos_sin(d['phi']) + machine_learning.angle_to_cos_sin(d['psi'])
              for d in feature_list]

    else:
      return [machine_learning.angle_to_cos_sin(d[0]) + machine_learning.angle_to_cos_sin(d[1])
              for d in feature_list]

  def learn(self, clf_type="OneClassSVM", transform_features=True):
    '''Learn the distribution with a machine learning classifier'''
    # Prepare the training data
  
    all_data = [(d['phi'], d['psi']) for d in self.feature_list]
    if transform_features:
      all_data = self.transform_features(all_data)
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
      self.clf = IsolationForest(max_samples=20000,
			contamination=0.01, random_state=np.random.RandomState(42))
      self.clf.fit(training_data)
   
    # Print Training results
    
    predictions = self.clf.predict(cv_data)
    print("{0}/{1} cross validation error.".format(len(predictions[-1 == predictions]), len(cv_data)))

    if clf_type == "OneClassSVM":
      print("{0} support vectors found.".format(len(self.clf.support_)))

  def predict(self, input_data, transform_features=True):
    '''Make a prediction for the input data with the machine learning classifier.
    input_data is a list of phi, psi angles.
    '''
    transformed_data = input_data

    if transform_features: 
        transformed_data = self.transform_features(transformed_data)
    
    return self.clf.predict(transformed_data)

  def calculate_space_reduction(self, transform_features=True):
    '''Calculate the space reduction power of the machine learning model.'''
    phis = np.random.uniform(-np.pi, np.pi, 10000)
    psis = np.random.uniform(-np.pi, np.pi, 10000)
    
    predictions = self.predict(list(zip(phis, psis)), transform_features=transform_features)
    print("The space is reduced by {0}.".format(len(predictions[1 == predictions]) / len(predictions)))

  def density_estimate(self, de_type="GaussianMixture", transform_features=True):
    '''Get a density estimation of the data.'''
    all_data = [(d['phi'], d['psi']) for d in self.feature_list]
    if transform_features:
      all_data = self.transform_features(all_data)
    n_data = len(all_data)
    training_data = all_data[0:int(0.7 * n_data)]
    test_data = all_data[int(0.7 * n_data):n_data]

    # Make some random data
    
    phis = np.random.uniform(-np.pi, np.pi, 10000)
    psis = np.random.uniform(-np.pi, np.pi, 10000)
    
    random_data = list(zip(phis, psis))
    if transform_features: 
        random_data = self.transform_features(random_data)

    if de_type == "GaussianMixture":
      self.de = mixture.BayesianGaussianMixture(n_components=100, covariance_type='full').fit(training_data)
      
      # Evalute the cumulative distribution functions of scores of test data

      test_scores = self.de.score_samples(test_data)
      values, base = np.histogram(test_scores, bins=40)
      cumulative = np.cumsum(values)

      for i in range(40):
        
        # Evaluate the space compression

        random_scores = self.de.score_samples(random_data)
        compress_coe = len(random_scores[random_scores > base[i]]) / len(random_scores)
          
        print('{0:.3f}\t{1}\t{2:.5f}\t{3:.5f}'.format(base[i], cumulative[i], cumulative[i] / len(test_data), compress_coe))


    elif de_type == "KernelDensity":
      params = {'bandwidth': np.logspace(-1, 1, 5)}
      grid = GridSearchCV(KernelDensity(), params)
      grid.fit(training_data) 
      self.de = grid.best_estimator_
