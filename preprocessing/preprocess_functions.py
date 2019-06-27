import pandas as pd
import os
import urllib
import xml.etree.ElementTree as ET

def partekSpecificPreprocessing(group_id, directory, classification):
    # INPUT: single raw partek file
    # OUTPUT:  -  cell_id_raw_name_map.csv (map of new cell_id with old cell_id)
    #          -  Gene_Symbol/* (every cell gene symbol list w/ expression value)
    #          -  all_expressed_genes.txt (list of all the probed genes)

    print("Starting partex specific preprocessing...")

    raw_file = '../raw_data/' + group_id + '.txt'

    if False and os.path.isdir(directory + group_id):
	    print("\"" + group_id + "\" directory already exists." \
                   + " Please delete existing directory to run partek specific preprocessing.")

    else:

        # os.mkdir(directory + group_id)

        print("Downloading raw data...")
        raw_data = pd.read_csv(raw_file,sep='\t')
        print("Successfully downloaded raw data.")

        directory = directory + group_id + "/"

        filtered_data = raw_data[raw_data['Classifications'] == classification]
        cell_id_raw_name_map = filtered_data[['Name']]
        cell_id_raw_name_map = cell_id_raw_name_map.rename(columns={'Name': 'raw_id'})
        cell_id_raw_name_map.to_csv(directory + \
                'cell_id_raw_name_map.csv', index_label='cell_id')

        print("Created cell_id to raw_name map.")

        long_filtered_data = filtered_data.transpose()
        cell_expression_array = long_filtered_data[8:]
       
        gene_list = cell_expression_array.index.tolist()

        
        with open(directory + 'all_expressed_genes.txt', 'w') as f:
            for gene in gene_list:
                f.write("%s\n" % gene)


        print("Created list of all expressed genes")

        os.mkdir(directory + 'Gene_Symbol')

        print("Making gene symbol files...")
        for cell_id in cell_expression_array.columns:
                print(cell_id)
                sc_data = cell_expression_array[[cell_id]]
                sc_data = sc_data.rename(columns={sc_data.columns[0] : 'expr_val'})
                sc_data.to_csv(directory+'Gene_Symbol/'
                        +group_id+'_'+str(cell_id)+'_gene_symbol.csv',index_label='gene_symbol')

        print("Made gene symbol files.")
        print("Partek specific preprocessing completed!")


def getEnsemblMap(group_id, directory):
    # INPUT: file containing all the gene symbols
    # OUTPUT: file containing all the gene symbols w/
    #         corresponding ensembl IDs

    print('Starting getEnsemblMap...')
    if os.path.exists(directory + group_id + '/gene_ensembl_map_all_probes.csv'):
	print("File gene_ensembl_map_all_probes.csv already exists." \
               + " Please delete existing file to get mapping.")
    else:
        
        directory = directory + group_id + '/'

        os.mkdir(directory + 'db2db_GeneToEns_Responses')
        
        # url query
        url_first = 'https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/' \
                    + 'biodbnetRestApi.xml?method=db2db&format=row&input=' \
                    + 'genesymbol&inputValues='
        url_last = '&outputs=ensemblproteinid&taxonId=10090'

        query_list = []

        chunk_size = 100

        url_genes = ''

        cur = 0

        f = open(directory + 'all_expressed_genes.txt','r')

        while True:
            
            gene = f.readline()
            url_genes += gene.replace('\n',',')
            cur += 1
            
            if (not gene) or cur%chunk_size == 0:
                
                if url_genes[-1:] == ',':
                    url_genes = url_genes[:-1]
                
                query_list.append(url_genes)
                
                url_genes = ''
                
                if not gene:
                    break
                
               
        print("Starting download from db2db.")
        print(str(len(query_list)) + " API calls.")
        count = 0
        for idx in range(len(query_list)):
            url_genes = query_list[idx]

            url = url_first + url_genes + url_last
            u = urllib.urlopen(url)
            response = u.read()

            with open(directory + 'db2db_GeneToEns_Responses/response_' + \
                      str(idx) + '.xml', 'w') as f:
                f.write(response)
                count+=1
                print(str(count) + "/" + str(len(query_list)) + " completed")

        response_list = os.listdir(directory + '/db2db_GeneToEns_Responses/')

        print("Beginning merge.")
        # initialize gene symbol to ensembl map
        ensembl_map = []
        
        for response_file in response_list:

            # create element tree object
            tree = ET.parse(directory + 'db2db_GeneToEns_Responses/' + response_file)

            # get root element
            root = tree.getroot()

            # fill map
            for item in root.findall('./'):
                for element in item:
                    if element.tag == 'InputValue':
                        gene_symbol = element.text
                    if element.tag == 'EnsemblProteinID':
                        if element.text == None:
                            best_ens_id = ""
                        else:
                            ens_id_list = element.text.split('//')
                            best_ens_id = ''
                            for val in range(len(ens_id_list)):
                                ens_id = ens_id_list[val]
                                if len(ens_id) > len(best_ens_id):
                                    best_ens_id = ens_id


                        ensembl_map.append((gene_symbol,best_ens_id))

        output_map_df = pd.DataFrame(ensembl_map,columns=['gene_symbol', 'ensembl_id'])
        output_map_df = output_map_df.iloc[output_map_df.gene_symbol.str.lower().argsort()]
        output_map_df.to_csv(directory + 'gene_ensembl_map_all_probes.csv',index=False)
        print("Completed merge.")

        print("getEnsemblMap completed!")


    return 0


def applyAllEnsemblMaps(group_id, directory):
    # INPUT: single gene symbol list file w/ expr val 
    # OUTPUT: single gene symbol and ensembl ID list w/ expr val

    print('Starting applyAllEnsemblMaps...')
    if os.path.isdir(directory + group_id + '/EnsemblID/'):
	print("\"" + group_id + "\"/EnsemblID/ directory already exists." \
               + " Please delete existing directory to run mapping.")

    else:

        directory = directory + group_id + '/'

        os.mkdir(directory + 'EnsemblID')

        gene_symbol_list = os.listdir(directory + 'Gene_Symbol')
        gene_ensembl_map = pd.read_csv(directory + 'gene_ensembl_map_all_probes.csv')
        gene_ensembl_dict = gene_ensembl_map.set_index('gene_symbol').to_dict()

        print("Create used subset of map in gene_ensembl_map.csv")
        # Computes the used gene ensembl map
        map_all_probes = pd.read_csv(directory + 'gene_ensembl_map_all_probes.csv')
        map_all_probes['ensembl_id'] = map_all_probes['gene_symbol'].map( \
                gene_ensembl_dict['ensembl_id'])

        map_all_probes = map_all_probes[map_all_probes['ensembl_id'].apply( \
            lambda x: not isinstance(x,float) and len(x)>1)]

        map_all_probes['ensembl_id'] = map_all_probes['ensembl_id'].apply( \
                lambda x: '10090.' + x)

        map_all_probes.to_csv(directory + 'gene_ensembl_map.csv', index=False)

        print("Finished creating gene_ensembl_map.csv")

        for gene_symbol_file in gene_symbol_list:

            gene_symbol_df = pd.read_csv(directory + 'Gene_Symbol/' + gene_symbol_file)
            gene_symbol_df['ensembl_id'] = gene_symbol_df['gene_symbol'].map( \
                                           gene_ensembl_dict['ensembl_id'])

            gene_symbol_df = gene_symbol_df[gene_symbol_df['ensembl_id'].apply( \
                                lambda x: not isinstance(x,float) and len(x)>1)]

            gene_symbol_df['ensembl_id'] = gene_symbol_df['ensembl_id'].apply( \
                                                        lambda x: '10090.' + x)

            output_filename = gene_symbol_file.replace('gene_symbol','ensembl')
            gene_symbol_df.to_csv(directory + 'EnsemblID/' + output_filename, index=False)

        print('applyAllEnsemblMaps completed!')

    return 0

 
def getNetworkData(group_id, directory):

    print('Starting getNetworkData...')
    if os.path.isdir('../network_data/' + group_id):

        print('network_data/' + group_id +'/ already exists. Please delete' \
               + ' existing directory to run mapping.')

    else: 

        os.mkdir('../network_data/' + group_id)

        ensembl_id_list = os.listdir(directory + group_id + '/EnsemblID')

        interactions = pd.read_csv(directory + group_id + '/ensembl_interactions.csv')
        interactions = interactions[['from', 'to']]

        
        for ensembl_id_file in ensembl_id_list:
            temp_interactions = interactions.copy()
            ens_id_df = pd.read_csv(directory + group_id + '/EnsemblID/' + ensembl_id_file)
            ens_id_df = ens_id_df[ens_id_df['expr_val'] > 0]
            ens_id_df = ens_id_df.drop(columns="gene_symbol")
            ens_id_df = ens_id_df[['ensembl_id','expr_val']]


            nodes = set(ens_id_df['ensembl_id'])
            in_from = temp_interactions['from'].apply(lambda x: x in nodes)
            in_to = temp_interactions['to'].apply(lambda x: x in nodes)
            temp_interactions = temp_interactions[in_from & in_to]
            edgelist_file = ensembl_id_file.replace('ensembl','edgelist')
            nodelist_file = ensembl_id_file.replace('ensembl','nodelist')
            temp_interactions.to_csv('../network_data/' + group_id + '/' + \
                    edgelist_file, index=False, header=False)
            ens_id_df.to_csv('../network_data/' + group_id + '/' + \
                    nodelist_file, index=False, header=False)

        print('getNetworkData completed!')


    return 0








