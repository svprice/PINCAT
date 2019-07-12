import os
import sys
import csv

import network_features




directory = '../network_data/giulia_t_lymphocytes/'

all_files = os.listdir(directory)

all_cell_ids = set(['_'.join(x.split('_')[:-1]) for x in all_files])

with open('giulia_t_lymphocytes_spatial_entropy.csv', 'w') as csvfile:

	fieldnames = ['nodelist','edgelist','num_nodes','num_edges', \
		    'giulia_spatial_entropy','runtime']

	writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
	writer.writeheader()

	for cell_id in all_cell_ids:
	    print(cell_id)

	    edgelist = directory + cell_id + '_edgelist.csv'
	    nodelist = directory + cell_id + '_nodelist.csv'

	    writer.writerow(network_features.giulia_spatial_entropy(edgelist,nodelist))


