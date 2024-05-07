import pandas as pd
import csv
import os
from rdkit import Chem

data_dir = './drug-files/'

drug_molecule_ids_dict = {}

MOLECULE_BOND_THRESHOLD = 19

for db_id in os.listdir(data_dir):
    molecule_file = os.path.join(data_dir, db_id)
    # print(molecule_file)
    if os.path.isfile(molecule_file):
        with Chem.SDMolSupplier(molecule_file) as suppl:
            molecule = [x for x in suppl if x is not None]
            drug_molecule_ids_dict[db_id[:-4]] = molecule[0]

common_mcs_file = 'common_mcs_antidepressants.txt'
common_mcs_dict = {}

with open(common_mcs_file, 'r') as fp:
    fp.readline()
    index = 0
    for mcs in fp:
        mcs_molecule = Chem.MolFromSmarts(mcs)
        bonds = mcs_molecule.GetNumBonds()
        if not common_mcs_dict.__contains__(mcs) and bonds > MOLECULE_BOND_THRESHOLD:
            common_mcs_dict[mcs] = index
            index = index + 1

feature_vector_length = len(common_mcs_dict)
print('MCS dictionary length: ', feature_vector_length)
print('Index: ', index)

dataset_header = []

for feature_number in range(feature_vector_length):
    dataset_header.insert(feature_number, 'feature_' + str(feature_number) + '_drug1')
    dataset_header.insert(feature_number + feature_vector_length, 'feature_' + str(feature_number) + '_drug2')

dataset_header.insert(2 * feature_vector_length + 1, 'Increase Activity Interaction')
dataset_header.insert(2 * feature_vector_length + 2, 'Decrease Activity Interaction')
dataset_header.insert(2 * feature_vector_length + 3, 'Increase Effect Interaction')
dataset_header.insert(2 * feature_vector_length + 4, 'Decrease Effect Interaction')
dataset_header.insert(2 * feature_vector_length + 5, 'Increase Efficacy Interaction')
dataset_header.insert(2 * feature_vector_length + 6, 'Decrease Efficacy Interaction')
dataset_header.insert(2 * feature_vector_length + 7, 'Other Interaction')

print(dataset_header)
print('Length of data header: ', len(dataset_header))

data = pd.DataFrame(columns=dataset_header)

feature_vector_list = []

source_file = open('drug_interaction_data_antidepressants.csv', encoding="utf8", newline='')
drug_interaction_file_reader = csv.reader(source_file, delimiter=",", quotechar='"')
n = 0

molecule_feature_dir = {}

common_mcs = common_mcs_dict.keys()
for drug_id, molecule_drug in drug_molecule_ids_dict.items():
    feature_vector = {}
    for mcs in common_mcs:
        mcs_index = common_mcs_dict.get(mcs)
        feature_column_name = dataset_header[mcs_index]
        mcs_molecule = Chem.MolFromSmarts(mcs)
        if molecule_drug.HasSubstructMatch(mcs_molecule):
            feature_vector[feature_column_name] = 1
        else:
            feature_vector[feature_column_name] = 0
    molecule_feature_dir[drug_id] = feature_vector

print('Molecule feature dir length: ', len(molecule_feature_dir))

for row_nr, row in enumerate(drug_interaction_file_reader):
    if n >= 20000:
        if n % 1000 == 0:
            print('Parsed %d interactions.' % n)
        if row_nr != 0:
            drug1 = row[0]
            drug2 = row[1]
            increase_activity_interaction = row[2]
            decrease_activity_interaction = row[3]
            increase_effect_interaction = row[4]
            decrease_effect_interaction = row[5]
            increase_efficacy_interaction = row[6]
            decrease_efficacy_interaction = row[7]
            other_interaction = row[8]
            # print("Drug 1: %s, Drug 2: %s, interaction: %s" % (drug1, drug2, interaction))
            # common_mcs = common_mcs_dict.keys()
            if molecule_feature_dir.__contains__(drug1) and molecule_feature_dir.__contains__(drug2):
                feature1 = molecule_feature_dir[drug1]
                feature2 = molecule_feature_dir[drug2]
                feature_vector = {**feature1, **feature2, 'Increase Activity Interaction': increase_activity_interaction,
                                  'Decrease Activity Interaction': decrease_activity_interaction,
                                  'Increase Effect Interaction': increase_effect_interaction,
                                  'Decrease Effect Interaction': decrease_effect_interaction,
                                  'Increase Efficacy Interaction': increase_efficacy_interaction,
                                  'Decrease Efficacy Interaction': decrease_efficacy_interaction,
                                  'Other Interaction': other_interaction}
                data = pd.concat([data, pd.DataFrame([feature_vector])], ignore_index=True)
            else:
                print('Interaction skipped cause sdf file not found: ', drug1, drug2)
    n += 1

train_test_file_name = "drug_interaction_train_test_ablation_20_antidepressants_2.csv"

drug_interaction_df = pd.read_csv(train_test_file_name)

data_frame_series = pd.concat([drug_interaction_df, data])

train_test_file_name_2 = "drug_interaction_train_test_ablation_20_antidepressants_final.csv"

print('Writing the first 20k interactions in file {}', train_test_file_name_2)
data_frame_series.to_csv(train_test_file_name_2, sep=',', encoding='utf-8')
print(data_frame_series.head(10))

# molecule_drug1 = drug_molecule_ids_dict.get(drug1)
# molecule_drug2 = drug_molecule_ids_dict.get(drug2)
# for mcs in common_mcs:
#     mcs_index = common_mcs_dict.get(mcs)
#     feature_column_name = dataset_header[mcs_index]
#     # print('Feature column drug1: ', feature_column_name)
#     feature2_column_name = dataset_header[mcs_index + feature_vector_length]
#     # print('Feature column drug2: ', feature2_column_name)
#     mcs_molecule = Chem.MolFromSmarts(mcs)
#     if molecule_drug1.HasSubstructMatch(mcs_molecule):
#         feature_vector[feature_column_name] = 1
#     else:
#         feature_vector[feature_column_name] = 0
#     if molecule_drug2.HasSubstructMatch(mcs_molecule):
#         feature_vector[feature2_column_name] = 1
#     else:
#         feature_vector[feature2_column_name] = 0