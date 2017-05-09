import numpy as np
import scipy
import scipy.spatial
import Bio.PDB as PDB


def angle(v1, v2):
    '''Return the angle between two vectors in radian.'''
    cos = max(-1, min(1, np.dot(normalize(v1), normalize(v2))))
    return np.arccos(cos)

def dihedral(p1, p2, p3, p4):
    '''Return the dihedral defined by 4 points in 
    range [-pi, pi].
    '''
    v1 = normalize(p2 - p1)
    v2 = normalize(p3 - p2)
    v3 = normalize(p4 - p3)

    n1 = normalize(np.cross(v1, v2))
    n2 = normalize(np.cross(v2, v3))

    c = np.dot(n1, n2)
    s = np.dot(v2, np.cross(n1, n2))

    return np.arctan2(s, c) 

def get_phi(chain, residue):
  '''Calculate the phi torsion of a residue.'''
  
  # Get the previous residue

  res_id = residue.get_id()
  prev_res = chain[res_id[1] - 1]
  
  prev_flag = prev_res.get_id()[0]
  if prev_flag == 'W' or prev_flag.startswith('H_'):
    raise Exception('Hetero residue type!')
 
  # Calculate the torsion

  c_prev = prev_res['C'].get_vector()
  n = residue['N'].get_vector()
  ca = residue['CA'].get_vector()
  c = residue['C'].get_vector()

  return PDB.calc_dihedral(c_prev, n, ca, c) 

def get_psi(chain, residue):
  '''Calculate the psi torsion of a residue.'''
  
  # Get the next residue

  res_id = residue.get_id()
  next_res = chain[res_id[1] + 1]
  
  next_flag = next_res.get_id()[0]
  if next_flag == 'W' or next_flag.startswith('H_'):
    raise Exception('Hetero residue type!')
  
  # Calculate the torsion

  n = residue['N'].get_vector()
  ca = residue['CA'].get_vector()
  c = residue['C'].get_vector()
  n_next = next_res['N'].get_vector()

  return PDB.calc_dihedral(n, ca, c, n_next)

def get_distance_matrix(atom_list):
  '''Get the distance matrix of a list of atoms.'''
  return scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(
      np.array([a.get_coord() for a in atom_list]), 'euclidean'))

def get_nearest_nonbonded_residues(model):
  '''Return a list of 2-tuples. Each the first element of each tuple is a
  residue and the second element is its nearest nonbonded residue.'''
 
  # Get all CA atoms which are used to be the center of residues

  ca_list = []
  
  for residue in model.get_residues():
    flag = residue.get_id()[0]
    if flag == 'W' or flag.startswith('H_'): continue
    
    for a in residue:
      if a.get_id() == 'CA':
        ca_list.append(a)

  ca_coords = [a.get_coord() for a in ca_list]

  # Make a KDTree for neighbor searching

  kd_tree = scipy.spatial.KDTree(ca_coords)

  # Find the nearest nonbonded neighbor of all residues

  nearest_nb_list = []

  for i in range(len(ca_list)):
    res1 = ca_list[i].get_parent()
    distance, indices = kd_tree.query(ca_coords[i], k=4)

    for j in range(1, 4):
      res2 = ca_list[indices[j]].get_parent() 

      if res1.get_parent().get_id() == res2.get_parent().get_id() \
         and (res1.get_id()[1] + 1 == res2.get_id()[1] \
         or res1.get_id()[1] - 1 == res2.get_id()[1]) : # Bonded residues
           continue
    
    nearest_nb_list.append((res1, res2))   

  return nearest_nb_list

def normalize(v):
  '''Normalize a numpy array.'''
  norm=np.linalg.norm(v)
  if norm==0: 
     return v
  return v/norm

def get_stub_matrix(p1, p2, p3):
  '''Get a matrix corresponding to a coordinate frame formed by 3 points.
     The origin is on p2, the y-axis is from p2 to p3; the z-axis is the
     cross product of vector p1-p2 and the y-axis; the x-axis is the product
     of y-axis and z-axis. The axis vectors will be the columns of the returned
     matrix.
  '''
  y = normalize(p3 - p2)
  z = normalize(np.cross(p1 - p2, y))
  x = np.cross(y, z)

  return np.matrix([x, y, z]).T

def get_residue_stub_matrix(residue):
  '''Constructure a coordinate frame on a residue. The origin is on the CA atom; 
  the y-axis is from CA to C; the z-axis is the cross product of the CA-N vector and the
  y-axis; the x-axis is thus defined by requiring the frame to be right handed.
  Return a 3x3 matrix and a vector that transforms coordinates in the local frame to 
  coordinates in the global frame.
  '''
  n = residue['N'].get_coord()
  ca = residue['CA'].get_coord()
  c = residue['C'].get_coord()

  return get_stub_matrix(n, ca, c), ca

def rotation_matrix_to_euler_angles(m):
  '''Return the euler angles corresponding to a rotation matrix.'''
  theta_x = np.arctan2(m[2][1], m[2][2])
  theta_y = np.arctan2(-m[2][0], np.sqrt(m[2][1]**2 + m[2][2]**2))
  theta_z = np.arctan2(m[1][0], m[0][0])

  return theta_x, theta_y, theta_z

def euler_angles_to_rotation_matrix(theta_x, theta_y, theta_z):
  '''Return the rotation matrix corresponding to 3 Euler angles.'''
  cx = np.cos(theta_x)
  sx = np.sin(theta_x)
  cy = np.cos(theta_y)
  sy = np.sin(theta_y)
  cz = np.cos(theta_z)
  sz = np.sin(theta_z)
  
  X = np.array([[  1,   0,   0],
                [  0,  cx, -sx],
                [  0,  sx,  cx]])
  Y = np.array([[ cy,   0,  sy],
                [  0,   1,   0],
                [-sy,   0,  cy]])
  Z = np.array([[ cz, -sz,   0],
                [ sz,  cz,   0],
                [  0,   0,   1]])

  return np.matmul(Z, np.matmul(Y, X))

def random_unit_vector(dim=3):
  '''Generate a random unit vector following the
  uniform distribution on the (dim - 1) dimension sphere.
  '''
  while True:
    v = np.random.normal(size=dim)
    if np.linalg.norm(v) > 0: break
  
  return normalize(v)

def random_rotation_matrix():
  '''Generate a random rotation matrix following the
  uniform distribution in SO(3).
  '''
  x = random_unit_vector()
  t = random_unit_vector()
  while np.linalg.norm(x - t) == 0:
    t = random_unit_vector()
  y = normalize(np.cross(x, t))
  z = np.cross(x, y)

  return np.array([x, y, z])

def random_euler_angles():
  '''Generate a random euler angles following the
  uniform distribution in SO(3).
  '''
  return rotation_matrix_to_euler_angles(random_rotation_matrix())

def rotation_matrix_from_axis_and_angle(u, theta):
  '''Calculate a rotation matrix from an axis and an angle.'''

  u = normalize(u)
  x = u[0]
  y = u[1]
  z = u[2]
  s = np.sin(theta)
  c = np.cos(theta)

  return np.array([[c + x**2 * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
                   [y * x * (1 - c) + z * s, c + y**2 * (1 - c), y * z * (1 - c) - x * s ],
                   [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z**2 * (1 - c) ]])

def cartesian_coord_from_internal_coord(p1, p2, p3, d, theta, tau):
  '''Calculate the cartesian coordinates of an atom from 
  three internal coordinates and three reference points.
  '''
  axis1 = np.cross(p1 - p2, p3 - p2)
  axis2 = p3 - p2

  M1 = rotation_matrix_from_axis_and_angle(axis1, theta - np.pi)
  M2 = rotation_matrix_from_axis_and_angle(axis2, tau)

  return p3 + d * np.dot(M2, np.dot(M1, normalize(p3 - p2)))

def cartesian_coord_to_internal_coord(p1, p2, p3, p):
  '''Cacluate internal coordinates of a point from its
  cartesian coordinates based on three reference points.
  '''
  return [np.linalg.norm(p - p3), angle(p - p3, p2 - p3), dihedral(p1, p2, p3, p)]

def calc_discrete_curvatures(p, neighbors):
  '''Calculate discrete curvatures of a point p
  given its neighbors in a cyclic order. Return the Gaussian
  curvature at the vertex and curvatures at the edges.
  '''
  n = len(neighbors)
  edges = [n - p for n in neighbors]
  
  # Cacluate the Gaussian curvature

  gaussian_curvature = 2 * np.pi

  for i in range(n):
    gaussian_curvature -= angle(edges[i], edges[(i + 1) % n])

  # Calculate edge curvatures

  edge_curvatures = []

  for i in range(n):
    edge_curvatures.append(angle(np.cross(edges[i], edges[(i + 1) % n]),
        np.cross(edges[(i + 1) % n], edges[(i + 2) % n])))

  return gaussian_curvature, edge_curvatures

