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

    with open(myod_group + '_spatial_entropy.csv', 'w') as csvfile:

        fieldnames = ['nodelist','edgelist','num_nodes','num_edges', \
                    'giulia_spatial_entropy','runtime']

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for cell_id in all_cell_ids:
            print(cell_id)

            edgelist = directory + myod_group + '_' + cell_id + '_edgelist.csv'
            nodelist = directory + myod_group + '_' + cell_id + '_nodelist.csv'

            writer.writerow(network_features.giulia_spatial_entropy(edgelist,nodelist))
        
    print('\n\n') 
