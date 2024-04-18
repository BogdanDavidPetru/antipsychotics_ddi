import pandas as pd
import csv
import os
from rdkit import Chem

data_dir = 'C:/Users/david.bogdan/master/disertatie/antipsychotic-drugbank/train_test'

drug_molecule_ids_dict = {}

for db_id in os.listdir(data_dir):
    molecule_file = os.path.join(data_dir, db_id)
    # print(molecule_file)
    if os.path.isfile(molecule_file):
        with Chem.SDMolSupplier(molecule_file) as suppl:
            molecule = [x for x in suppl if x is not None]
            drug_molecule_ids_dict[db_id[:-4]] = molecule[0]

common_mcs_file = 'common_mcs_antipsychotics.txt'
common_mcs_dict = {}

with open(common_mcs_file, 'r') as fp:
    fp.readline()
    index = 0
    for mcs in fp:
        if not common_mcs_dict.__contains__(mcs):
            common_mcs_dict[mcs] = index
            index = index + 1

feature_vector_length = len(common_mcs_dict)
print('MCS dictionary length: ', feature_vector_length)
print('Index: ', index)

dataset_header = []

for feature_number in range(feature_vector_length):
    dataset_header.insert(feature_number, 'feature_' + str(feature_number) + '_drug1')
    dataset_header.insert(feature_number + feature_vector_length, 'feature_' + str(feature_number) + '_drug2')

dataset_header.insert(2 * feature_vector_length + 1, 'interaction')

print(dataset_header)
print('Length of data header: ', len(dataset_header))

data = pd.DataFrame(columns=dataset_header)

feature_vector_list = []

source_file = open('drug_interaction_data_not_valid_antipsychotics.csv', encoding="utf8", newline='')
drug_interaction_file_reader = csv.reader(source_file, delimiter=",", quotechar='"')

for row_nr, row in enumerate(drug_interaction_file_reader):
    if row_nr != 0:
        drug1 = row[0]
        drug2 = row[1]
        interaction = row[2]
        # print("Drug 1: %s, Drug 2: %s, interaction: %s" % (drug1, drug2, interaction))
        common_mcs = common_mcs_dict.keys()
        feature_vector = {}
        molecule_drug1 = drug_molecule_ids_dict.get(drug1)
        molecule_drug2 = drug_molecule_ids_dict.get(drug2)
        for mcs in common_mcs:
            mcs_index = common_mcs_dict.get(mcs)
            feature_column_name = dataset_header[mcs_index]
            # print('Feature column drug1: ', feature_column_name)
            feature2_column_name = dataset_header[mcs_index + feature_vector_length]
            # print('Feature column drug2: ', feature2_column_name)
            mcs_molecule = Chem.MolFromSmarts(mcs)
            if molecule_drug1.HasSubstructMatch(mcs_molecule):
                feature_vector[feature_column_name] = 1
            else:
                feature_vector[feature_column_name] = 0
            if molecule_drug2.HasSubstructMatch(mcs_molecule):
                feature_vector[feature2_column_name] = 1
            else:
                feature_vector[feature2_column_name] = 0
        feature_vector['interaction'] = interaction
        data = pd.concat([data, pd.DataFrame([feature_vector])], ignore_index=True)

train_test_file_name = "drug_interaction_train_test_data_not_valid_antipsychotics.csv"
data.to_csv(train_test_file_name, sep=',', encoding='utf-8')
print(data.head(10))
