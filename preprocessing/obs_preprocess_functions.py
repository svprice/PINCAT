''' OBSOLETE
def downloadEnsemblMapData(group_id, directory):
    # INPUT: file containing all the gene symbols
    # OUTPUT: file containing all the gene symbols w/
    #         corresponding ensembl IDs

    print('Starting downloadEnsemblMapData...')
    if os.path.exists(directory + group_id + '/db2db_GeneToEns_Responses'):
	print("Directory db2db_GeneToEns_Responses already exists." \
               + " Please delete existing directory to download data.")
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

    print('downloadEnsemblMapData completed!')

def mergeToCreateEnsemblMap(group_id, directory):

    print("Starting mergeToCreateEnsemblMap...")

    if os.path.exists(directory + group_id + '/gene_ensembl_map.csv'):
	print("File gene_ensembl_map.csv already exists." \
               + " Please delete existing file to merge data.")
    else:

        directory = directory + group_id + '/'

        response_list = os.listdir(directory + '/db2db_GeneToEns_Responses/')

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
                            ensembl_map.append((gene_symbol,""))
                        else:
                            ens_id_list = element.text.split('//')
                            for ens_id in ens_id_list:
                                ensembl_map.append((gene_symbol,ens_id))



        output_map_df = pd.DataFrame(ensembl_map,columns=['gene_symbol', 'ensembl_id'])
        output_map_df = output_map_df.iloc[output_map_df.gene_symbol.str.lower().argsort()]

        output_map_df['ensembl_id'] = output_map_df['ensembl_id'].apply( \
                lambda x: '10090.' + x)

        output_map_df.to_csv(directory + 'gene_ensembl_map.csv',index=False)

        print("mergeToCreateEnsemblMap completed!")

    return 0
'''

