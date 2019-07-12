import os
import sys
import csv

import network_features

myod_folders = ['partek_myod_old', 'partek_myod_young', 'partek_myod_postnatal']


for myod_group in myod_folders:

    print(myod_group)

    directory = '../network_data/' + myod_group + '/'

    all_files = os.listdir(directory)

    all_cell_ids = set([x.split('_')[-2] for x in all_files])

    with open(myod_group + '_summary_analysis.csv', 'w') as csvfile:

        fieldnames = ['cell_id',
                      'num_nodes',\
                      'num_edges', \
                      'average_node_degree', \
                      'num_connected_components', \
                      'node_connectivity', \
                      'average_clustering', \
                      'max_clique_size', \
                      'max_independent_set_size', \
                      'min_vertex_cover_size', \
                      'degree_assortativity_coefficient', \
                      'average_neighbor_degree', \
                      'diameter', \
                      'wiener_index', \
                      'runtime', \
                      'edgelist']
                      

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for cell_id in all_cell_ids:
            print(cell_id)

            edgelist = directory + myod_group + '_' + cell_id + '_edgelist.csv'
            nodelist = directory + myod_group + '_' + cell_id + '_nodelist.csv'

            output = network_features.get_summary(edgelist)
            output['cell_id'] = cell_id
            writer.writerow(output)
        
    print('\n\n') 
