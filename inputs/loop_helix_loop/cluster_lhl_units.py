#!/usr/bin/env python3
'''Cluster the LHL units 
Dump the lhl clusters to a clusters directory 

Usage:
    ./cluster_lhl_units.py lhl_info_path
'''

import os
import sys
import json

import numpy as np


def lhl_units_in_same_cluster(lhl_unit1, lhl_unit2, ca_cutoff):
    '''Return true if the two LHL units are in a same cluster.'''
    if lhl_unit1['pdb_file'] == lhl_unit2['pdb_file']:
        return False

    if np.linalg.norm(lhl_unit1['pre_anchor_ca'] - lhl_unit2['pre_anchor_ca']) > ca_cutoff:
        return False

    if np.linalg.norm(lhl_unit1['post_anchor_ca'] - lhl_unit2['post_anchor_ca']) > ca_cutoff:
        return False

    return True

def get_cluster_for_lhl_unit_recursive(lhl_unit_id, lhl_edges, lhl_units_clustered_mask, cluster_list):
    '''Find the cluster for the lhl unit.
    '''
    cluster_list.append(lhl_unit_id)

    lhl_units_clustered_mask[lhl_unit_id] = True

    for i in lhl_edges[lhl_unit_id]:
        if lhl_units_clustered_mask[i]:
            continue

        get_cluster_for_lhl_unit_recursive(i, lhl_edges, lhl_units_clustered_mask, cluster_list)

def cluster_lhl_units(lhl_info_path, ca_cutoff=2):
    '''Cluster the LHL units.
    Two LHL units are regarded as in the same cluster if
    their pre_anchor_cas and post_anchor_cas are within a
    cut off and they are from different pdb files.
    '''
    lhl_infos = []

    # Load the lhl units

    for lhl_info_file in os.listdir(lhl_info_path):

        with open(os.path.join(lhl_info_path, lhl_info_file), 'r') as f:
            lhl_info = json.load(f)

            lhl_infos += lhl_info

    # Convert the coordinates to numpy arrays 

    for i in range(len(lhl_infos)):
        lhl_infos[i]['pre_anchor_ca'] = np.array(lhl_infos[i]['pre_anchor_ca'])
        lhl_infos[i]['post_anchor_ca'] = np.array(lhl_infos[i]['post_anchor_ca'])

    # Find edges between the LHL units

    lhl_edges = []

    for i in range(len(lhl_infos)):
        lhl_edges.append([])

        for j in range(i + 1, len(lhl_infos)):

            if lhl_units_in_same_cluster(lhl_infos[i], lhl_infos[j], ca_cutoff):
                lhl_edges[i].append(j)

    # Get the clusters of lhl units 

    lhl_units_clustered_mask = [False] * len(lhl_infos)

    clusters = []

    for i in range(len(lhl_infos)):
        if lhl_units_clustered_mask[i]:
            continue

        clusters.append([])
        get_cluster_for_lhl_unit_recursive(i, lhl_edges, lhl_units_clustered_mask, clusters[-1])

    # Reload the lhl units

    lhl_infos = []

    for lhl_info_file in os.listdir(lhl_info_path):

        with open(os.path.join(lhl_info_path, lhl_info_file), 'r') as f:
            lhl_info = json.load(f)

            lhl_infos += lhl_info

    # Dump the clusters
    
    os.makedirs('clusters', exist_ok=True)
    clusters_sorted = sorted(clusters, key=lambda c:len(c), reverse=True)

    for i in range(len(clusters_sorted)):
        c = [lhl_infos[j] for j in clusters_sorted[i]]
       
        with open(os.path.join('clusters', '{0}.json'.format(i)), 'w') as f:
            json.dump(c, f)           

if __name__ == '__main__':
    lhl_info_path = sys.argv[1]

    cluster_lhl_units(lhl_info_path)
