import os
import sys
import csv

import network_features

directory = '../network_data/subpopulations_in_progress/'

all_dirs = os.listdir(directory)

dir_cur = 1
total_dirs = str(len(all_dirs))

for class_dir in all_dirs:
    
    all_files = os.listdir(directory + class_dir)

    all_cell_ids = set(['_'.join(x.split('_')[:-1]) for x in all_files])
    all_cell_ids.remove('cell_id_raw_name')

    print('\n\nThis is directory ' + str(dir_cur) + ' out of ' + total_dirs)
    dir_cur += 1 
    print('SUBPOPULATION: ' + class_dir)

    total = str(len(all_cell_ids))
    print(total + ' to compute\n')

    with open('subpopulation_analysis/subpop_' + class_dir + '.csv', 'w') as csvfile:

        fieldnames = ['subpopulation','cell_num', \
                      'num_nodes','num_edges', \
                      'giulia_spatial_entropy','runtime', \
                      'nodelist', 'edgelist']

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        cur = 0

        for cell_id in all_cell_ids:
                
            print(str(cur) + '/' + total + ' completed')
                
            cur+=1
            
            edgelist = directory + class_dir + '/' + cell_id + '_edgelist.csv'
            nodelist = directory + class_dir + '/' + cell_id + '_nodelist.csv'

            output = network_features.giulia_spatial_entropy(edgelist,nodelist)

            subpop = '_'.join(cell_id.split('_')[:-1])
            cell_num = cell_id.split('_')[-1]
                    
            output['subpopulation'] = subpop
            output['cell_num'] = cell_num
                    
            writer.writerow(output)

    print('Completed ' + class_dir)
