import pandas as pd
import os
import urllib
import xml.etree.ElementTree as ET

def partekSpecificPreprocessing(raw_file, output_directory, classification):
    # INPUT: single raw partek file
    # OUTPUT:  -  cell_id_raw_name_map.csv (map of new cell_id with old cell_id)
    #          -  Gene_Symbol/* (every cell gene symbol list w/ expression value)
    #          -  all_expressed_genes.txt (list of all the probed genes)

    raw_data = pd.read_csv(raw_file,sep='\t')

    filtered_data = raw_data[raw_data['Classifications'] == classification]
    cell_id_raw_name_map = filtered_data[['Name']]
    cell_id_raw_name_map = cell_id_raw_name_map.rename(columns={'Name': 'raw_id'})
    cell_id_raw_name_map.to_csv(output_directory + 'cell_id_raw_name_map.csv')

    long_filtered_data = filtered_data.transpose()
    cell_expression_array = long_filtered_data[8:]
   
    gene_list = cell_expression_array.index.tolist()

    with open(output_directory + 'all_expressed_genes.txt', 'w') as f:
        for gene in gene_list:
            f.write("%s\n" % gene)


    group_id = raw_file.replace('.csv','')

    os.mkdir(output_directory + 'Gene_Symbol')

    for cell_id in cell_expression_array.columns:
            sc_data = cell_expression_array[[cell_id]]
            sc_data = sc_data.rename(columns={sc_data.columns[0] : 'expr_val'})
            sc_data.to_csv(output_directory+'Gene_Symbol/'
                    +group_id+'_'+str(cell_id)+'_gene_symbol.csv',index_label='gene_symbol')


    # will output print data to bash script, but not return data
     


def madeUpFunction(a,b,c):
    print(b+a+c)
    print("this is the returned")


def getEnsemblMap(all_gene_symbol_file):
    # INPUT: file containing all the gene symbols
    # OUTPUT: file containing all the gene symbols w/
    #         corresponding ensembl IDs

    return 0


def applyEnsemblMap(gene_symbol_file, ensembl_map):
    # INPUT: single gene symbol list file w/ expr val 
    # OUTPUT: single gene symbol and ensembl ID list w/ expr val

    return 0

 


