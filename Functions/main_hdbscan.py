import argparse
import numpy as np
import json
import os
import scipy
import hdbscan

parser = argparse.ArgumentParser() 
parser.add_argument('json_filename', help='the name of json file')  
args = parser.parse_args()  

if args.json_filename:
    settings = json.load(open(args.json_filename))

data_folder = settings['data_folder']

distance_matrix = np.load(os.path.join(data_folder, 'DistanceMatrix.npy'))

mdl = hdbscan.HDBSCAN(
    min_cluster_size=settings['min_cluster_size'],
    metric=settings['metric'],
    max_cluster_size=settings['max_cluster_size'],
    min_samples=settings['min_samples'],
    cluster_selection_epsilon=settings['cluster_selection_epsilon']
)

idx = mdl.fit_predict(distance_matrix)
np.save(os.path.join(data_folder, 'ClusterIndices.npy'), idx)

Z = mdl.single_linkage_tree_.to_numpy()
Z_matlab = scipy.cluster.hierarchy.to_mlab_linkage(Z)
np.save(os.path.join(data_folder, 'LinkageMatrix.npy'), Z_matlab)

