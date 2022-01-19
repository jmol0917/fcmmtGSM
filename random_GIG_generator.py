"""
This script requires a .csv file with cofactor Boolean data and clear synthesis and requirements headings, (path_in) and a list of strain names (agora_files) to generate random
GIGs of a number of members equals to the value of number_of_members.

The number of random GIGs this script generates is equals to the value specified in the run_limit variable
"""

import pandas as pd
import os
import random
from os import listdir
from os.path import isfile, join

path_in = ''
agora_files = ''
cofactor_data = pd.read_csv(path_in, index_col=0)

# Creates a list of bacteria names (models) located in the path_in directory
models_in = [f for f in listdir(agora_files) if isfile(join(agora_files, f))]
models_in = [os.path.splitext(f)[0] for f in models_in]

run_limit = 1000
run = 0
number_of_strains = 9

while run < run_limit:
    i = 0
    
    test_gig = {
    'GIG': [],
    'Node': [],
    'Strains': []
    }
    
    while i < number_of_strains:
        r_number = random.randint(0, 815)
        test_gig['GIG'].append(run)
        test_gig['Node'].append(i)
        test_gig['Strains'].append(models_in[r_number])
        i += 1
    
    network_file = pd.DataFrame(test_gig)
    network_file = network_file.set_index('GIG')
    print(network_file)
    
    def cartesian_product(*lists):
        if len(lists) == 1:
            return [(x0,) for x0 in lists[0]]
        else:
            return [(x0,) + t1 for x0 in lists[0] for t1 in cartesian_product(*lists[1:])]

    list_of_GIGs = []


    for index, row in network_file.iterrows():
        if index not in list_of_GIGs:
            list_of_GIGs.append(index)

    for GIG in list_of_GIGs:
        version = 0
        strains_by_node = {}

        for index, row in network_file.iterrows():
            if index == GIG:

                list_of_strains_by_species = []
                strain = row['Strains']
                node = row['Node']
                if node not in strains_by_node:
                    strains_by_node[node] = []

                if ' ' or '-' in strain:
                    strain = strain.replace(' ', '_')
                    strain = strain.replace('-', '_')

                for model in models_in:
                    vit_requirements_list = []
                    if strain in model:
                        strains_by_node[node].append([node, model])


        lists_of_species_by_node = []

        for node in strains_by_node:
            pairs = strains_by_node[node]
            if len(pairs) == 1:
                for node2 in strains_by_node:
                    if node != node2:
                        pairs2 = strains_by_node[node2]
                        if pairs[0][1] == pairs2[0][1]:
                            print('Warning: two nodes contain the same single strain')
                            print(node, node2)
                            print('Modify input document to obtain accurate results.')

        for node in strains_by_node:
            if len(strains_by_node[node]) > 0:
                lists_of_species_by_node.append(strains_by_node[node])
            else:
                print(node, 'has no corresponding strains in the AGORA set')
        
        cartesian_product_of_strains = cartesian_product(*lists_of_species_by_node)

        version = 0
        for combination in cartesian_product_of_strains:

            cofactor_req_syn = {
                'Strain': [],
                'Node':[],
                'B1': [],
                'B2': [],
                'B3': [],
                'B5': [],
                'B6': [],
                'B9': [],
                'K': []
                }

    #         print(combination)
            strains_in_combination = []
            for pair in combination:
                if pair[1] not in strains_in_combination:
                    strains_in_combination.append(pair[1])


            if len(combination) == len(strains_in_combination):
    #             print(combination)
                for strain in combination:
    #                 print(strain)
                    vit_requirements_list = []
                    for index2, row2 in cofactor_data.iterrows():
                        if index2 == strain[1]:
                            cofactor_req_syn['Strain'].append(strain[1])
                            cofactor_req_syn['Node'].append(strain[0])
                            if int(row2.loc['B1_req']) > 0:
                                cofactor_req_syn['B1'].append('R')
                            elif int(row2.loc['B1_syn']) > 0:
                                cofactor_req_syn['B1'].append('S')
                            else:
                                cofactor_req_syn['B1'].append('-')
                            if int(row2.loc['B2_req']) > 0:
                                cofactor_req_syn['B2'].append('R')
                            elif int(row2.loc['B2_syn']) > 0:
                                cofactor_req_syn['B2'].append('S')
                            else:
                                cofactor_req_syn['B2'].append('-')
                            if int(row2.loc['B3_req']) > 0:
                                cofactor_req_syn['B3'].append('R')
                            elif int(row2.loc['B3_syn']) > 0:
                                cofactor_req_syn['B3'].append('S')
                            else:
                                cofactor_req_syn['B3'].append('-')
                            if int(row2.loc['B5_req']) > 0:
                                cofactor_req_syn['B5'].append('R')
                            elif int(row2.loc['B5_syn']) > 0:
                                cofactor_req_syn['B5'].append('S')
                            else:
                                cofactor_req_syn['B5'].append('-')
                            if int(row2.loc['B6_req']) > 0:
                                cofactor_req_syn['B6'].append('R')
                            elif int(row2.loc['B6_syn']) > 0:
                                cofactor_req_syn['B6'].append('S')
                            else:
                                cofactor_req_syn['B6'].append('-')
                            if int(row2.loc['B9_req']) > 0:
                                cofactor_req_syn['B9'].append('R')
                            elif int(row2.loc['B9_syn']) > 0:
                                cofactor_req_syn['B9'].append('S')
                            else:
                                cofactor_req_syn['B9'].append('-')
                            if int(row2.loc['K_req']) > 0:
                                cofactor_req_syn['K'].append('R')
                            elif int(row2.loc['K_syn']) > 0:
                                cofactor_req_syn['K'].append('S')
                            else:
                                cofactor_req_syn['K'].append('-')

                cofactor_req_syn_df = pd.DataFrame(cofactor_req_syn)
                cofactor_req_syn_df = cofactor_req_syn_df.set_index('Strain')
                cofactor_req_syn_df.to_csv('C:/Users/jpmo_/Dropbox (Sydney Uni)/Functional characterisation of the gut microbiome/test_networks/1000/GIG' + str(GIG) + '_version' + str(version) + '.csv')
                version += 1
    run += 1
