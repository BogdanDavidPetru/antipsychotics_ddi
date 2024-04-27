import os
import csv
import pandas as pd


class Triplets:
    def __init__(self):
        self.list = []

    def __len__(self):
        return len(self.list)

    def add(self, sujet, predicat, objet):
        if {"subject": sujet, "predicat": predicat, "object": objet} not in self.list:
            self.list.append({"subject": sujet,
                              "predicate": predicat,
                              "object": objet})


# not_researched_drugs = {12867, 12710, 12401, 13256, 13273, 13552, 13665, 13784, 13791, 13841, 13523, 4872, 4888, 8927, 11736, 12093, 12273, 12518, 12958, 13213, 13382, 13403, 13554, 13557, 13676, 372, 11540, 13420, 6077, 9266, 14185, 14651, 16021, 5766, 17056, 12543}

drug_name_file = './drug-files/info/drug_agent_list_drugbank.txt'

drug_names = set()
with open(drug_name_file, 'r') as fp:
    for drug in fp:
        drug_names.add(drug.strip())

interactions = []

known_ddi_file = ('C:/Users/david.bogdan/master/disertatie/oregano/oregano-master/oregano-master/Integration'
                  '/Integration V2.1/DrugBank/interaction_drugs_drugbank.tsv')

DDI = pd.read_csv(known_ddi_file, sep="\t", engine="python", names=["subject", "predicate", "object"])
# left_ddi_info = {}
# right_ddi_info = {}
DRUGS = {}
DDI_list = Triplets()

n = 0
identified_ddi = set()
for index in range(len(DDI["subject"])):
    if n % 100000 == 0:
        print(n, " / ", len(DDI["subject"]))
    if DDI["predicate"][index] == "is":
        if drug_names.__contains__(DDI["object"][index]):
            DRUGS[DDI["object"][index]] = DDI["subject"][index]
    else:
        if drug_names.__contains__(DDI["predicate"][index]) and drug_names.__contains__(DDI["subject"][index]):
            drug1 = DDI["subject"][index]
            drug2 = DDI["predicate"][index]
            DDI_list.add(drug1, drug2, DDI["object"][index])
            if not identified_ddi.__contains__((drug2, drug)):
                identified_ddi.add((drug1, drug2))
    n += 1

print("Parsing drug-drug interactions")

n = 0
increase_activity_interaction = set()
decrease_activity_interaction = set()

increase_effect_interaction = set()
decrease_effect_interaction = set()

increase_efficacy_interaction = set()
decrease_efficacy_interaction = set()

another_interaction = set()

filename = "drug_interaction_data_antidepressants.csv"

for ddi in DDI_list.list:
    if n % 10000 == 0:
        print(n, " / ", len(DDI_list.list))
    drug1 = DRUGS[ddi["predicate"]]
    drug2 = DRUGS[ddi["subject"]]
    if "may" in ddi["object"]:
        if "increase" in ddi["object"]:
            activity = ddi["object"].split("increase")[1].split("of")[0]
            if not increase_activity_interaction.__contains__((drug2, drug1)):
                increase_activity_interaction.add((drug1, drug2))
        else:
            activity = ddi["object"].split("decrease")[1].split("of")[0]
            if not decrease_activity_interaction.__contains__((drug2, drug1)):
                decrease_activity_interaction.add((drug1, drug2))
    elif "The risk or severity of" in ddi["object"]:
        effect = ddi["object"].split("of")[1].split("can be")[0]
        if "increased" in ddi["object"]:
            if not increase_effect_interaction.__contains__((drug2, drug1)):
                increase_effect_interaction.add((drug1, drug2))
        else:
            if not decrease_effect_interaction.__contains__((drug2, drug1)):
                decrease_effect_interaction.add((drug1, drug2))
    elif "The therapeutic efficacy of" in ddi["object"]:
        if "increased" in ddi["object"]:
            if not increase_efficacy_interaction.__contains__((drug2, drug1)):
                increase_efficacy_interaction.add((drug1, drug2))
        else:
            if not decrease_efficacy_interaction.__contains__((drug2, drug1)):
                decrease_efficacy_interaction.add((drug1, drug2))
    else:
        if not another_interaction.__contains__((drug2, drug1)):
            another_interaction.add((drug1, drug2))
    n += 1

drugs_without_interaction = set()

header = ['Drug1', 'Drug2', 'Increase Activity Interaction', 'Decrease Activity Interaction', 'Increase Effect '
                                                                                              'Interaction',
          'Decrease Effect Interaction', 'Increase Efficacy Interaction', 'Decrease Efficacy Interaction',
          'Other Interaction']
with open(filename, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(header)
    for index, drug in enumerate(drug_names):
        for index2, drug2 in enumerate(drug_names):
            if index2 > index:
                increase_activity_interaction_var = 0
                decrease_activity_interaction_var = 0
                increase_effect_interaction_var = 0
                decrease_effect_interaction_var = 0
                increase_efficacy_interaction_var = 0
                decrease_efficacy_interaction_var = 0
                other_interaction_var = 0
                if DRUGS.__contains__(drug) and DRUGS.__contains__(drug2):
                    drug_id = DRUGS[drug]
                    drug2_id = DRUGS[drug2]
                    if increase_activity_interaction.__contains__((drug_id, drug2_id)) or increase_activity_interaction.__contains__((drug2_id, drug_id)):
                        increase_activity_interaction_var = 1
                    if decrease_activity_interaction.__contains__((drug_id, drug2_id)) or decrease_activity_interaction.__contains__((drug2_id, drug_id)):
                        decrease_activity_interaction_var = 1
                    if increase_effect_interaction.__contains__((drug_id, drug2_id)) or increase_effect_interaction.__contains__((drug2_id, drug_id)):
                        increase_effect_interaction_var = 1
                    if decrease_effect_interaction.__contains__((drug_id, drug2_id)) or decrease_effect_interaction.__contains__((drug2_id, drug_id)):
                        decrease_effect_interaction_var = 1
                    if increase_efficacy_interaction.__contains__((drug_id, drug2_id)) or increase_efficacy_interaction.__contains__((drug2_id, drug_id)):
                        increase_efficacy_interaction_var = 1
                    if decrease_efficacy_interaction.__contains__((drug_id, drug2_id)) or decrease_efficacy_interaction.__contains__((drug2_id, drug_id)):
                        decrease_efficacy_interaction_var = 1
                    if another_interaction.__contains__((drug_id, drug2_id)) or another_interaction.__contains__((drug2_id, drug_id)):
                        other_interaction_var = 1
                    csvwriter.writerow([drug_id, drug2_id, increase_activity_interaction_var, decrease_activity_interaction_var, increase_effect_interaction_var, decrease_effect_interaction_var, increase_efficacy_interaction_var, decrease_efficacy_interaction_var, other_interaction_var])
                else:
                    if not DRUGS.__contains__(drug):
                        drugs_without_interaction.add(drug)
                    if not DRUGS.__contains__(drug2):
                        drugs_without_interaction.add(drug2)

print('Number of drugs taken into consideration for dataset: ', len(drug_names))
#
print("Interaction with increased activity: ", len(increase_activity_interaction))
print("Interaction with decreased activity: ", len(decrease_activity_interaction))

print("Interaction with increased effect: ", len(increase_effect_interaction))
print("Interaction with decreased effect: ", len(decrease_effect_interaction))

print("Interaction with increased efficacy: ", len(increase_efficacy_interaction))
print("Interaction with decreased efficacy: ", len(decrease_efficacy_interaction))

print("Drugs without interaction: ", drugs_without_interaction)

# print("Other interactions: ", another_interaction)
# print("--------------------DRUGS------------------")
# print(len(DRUGS))
# print(DRUGS)

# identified_ddi_file_name = "antidepressants_ddi.txt"

# identified_ddi_file = open(identified_ddi_file_name, "a")


# print("--------------------DDI Set------------------")
#
# for ddi in identified_ddi:
#     identified_ddi_file.write(str(ddi))
# identified_ddi_file.close()
# print("--------------------DDI List------------------")

# print(len(DDI_list))
# print(DDI_list.list)
# drugs_with_interactions = set()
# with open(known_DDI_file, 'r') as fp:
#     fp.readline()
#     for line in fp:
#         # print(line)
#         sptlist = line.split()  # line.strip().split('/t')
#
#         left_drug = sptlist[0].strip()
#         right_drug = sptlist[1].strip()
#         interaction_type = sptlist[2].strip()
#         if drug_ids.__contains__(left_drug) and drug_ids.__contains__(right_drug):
#             interactions.append((left_drug, right_drug, interaction_type))
#             drugs_with_interactions.add(left_drug)
#             drugs_with_interactions.add(right_drug)
#
# for interaction in interactions:
#     print("Drug 1: %s, Drug 2: %s, interaction: %s" % (interaction[0], interaction[1], interaction[2]))
#
# print("Interactions Size: ", len(interactions))
#
# filename = "drug_interaction_data_not_valid_antipsychotics.csv"
# header = ['Drug1', 'Drug2', 'Interaction']
#
# already_parsed = set()
#
# with open(filename, 'w', newline='') as csvfile:
#     csvwriter = csv.writer(csvfile)
#     csvwriter.writerow(header)
#     for index, drug in enumerate(drug_ids):
#         for index2, drug2 in enumerate(drug_ids):
#             if index2 > index:
#                 interaction_found = False
#                 for (left_drug, right_drug, interaction_type) in interactions:
#                     if (left_drug == drug and right_drug == drug2) or (left_drug == drug2 and right_drug == drug):
#                         interaction_found = True
#                         break
#                 if not already_parsed.__contains__((drug, drug2)) and not already_parsed.__contains__((drug2, drug)):
#                     if interaction_found:
#                         csvwriter.writerow([drug, drug2, 1])
#                     else:
#                         csvwriter.writerow([drug, drug2, 0])
#                     already_parsed.add((drug, drug2))
#
# print('Number of drugs taken into consideration for dataset: ', len(drug_ids))
# drugs_without_interactions = set()
#
# for drug_id in drug_ids:
#     if not drugs_with_interactions.__contains__(drug_id):
#         drugs_without_interactions.add(drug_id)
#
# print("Drugs without known interaction: ", len(drugs_without_interactions))
# print(drugs_without_interactions)
# print(drug_ids)
