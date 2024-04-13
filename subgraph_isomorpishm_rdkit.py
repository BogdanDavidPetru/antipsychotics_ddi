import os
import time
from rdkit import Chem
from rdkit.Chem import rdFMCS

data_dir = 'C:/Users/david.bogdan/master/disertatie/antipsychotic-drugbank/all'

drug_molecule = []

for mol in os.listdir(data_dir):
    full_path = os.path.join(data_dir, mol)
    if os.path.isfile(full_path):
        with Chem.SDMolSupplier(full_path) as suppl:
            ms = [x for x in suppl if x is not None]
            drug_molecule.append(ms)

print(len(drug_molecule))

mcs_set = []
start = time.time()

for index, molecule in enumerate(drug_molecule):
    for index2, molecule2 in enumerate(drug_molecule):
        if index2 > index:
            for mol in molecule:
                for mol2 in molecule2:
                    mcs = rdFMCS.FindMCS([mol, mol2])
                    db_id_1 = mol.GetProp('DRUGBANK_ID')
                    db_id_2 = mol2.GetProp('DRUGBANK_ID')
                    if mcs.numAtoms == 0 and mcs.numBonds == 0:
                        print("The result is an empty subgraph for (%s, %s)", db_id_1, db_id_2)
                    else:
                        mcs_set.append(mcs)
                    if len(mcs.smartsString) == 0:
                        print("The smartString is empty for (%s, %s)", db_id_1, db_id_2)

print(len(mcs_set))
f = open("./common_mcs.txt", "a")

for mcs in mcs_set:
    f.write(mcs.smartsString + "\n")
f.close()
# with Chem.SDMolSupplier('C:/Users/david.bogdan/master/disertatie/antipsychotic-drugbank/DB00372.sdf') as suppl:
#     ms = [x for x in suppl if x is not None]
#
# with Chem.SDMolSupplier('C:/Users/david.bogdan/master/disertatie/antipsychotic-drugbank//DB00246.sdf') as suppl2:
#     ms2 = [x for x in suppl2 if x is not None]
#
# start = time.time()
# for molecule in ms:
#     for second_molecule in ms2:
#
#         mcs = rdFMCS.FindMCS([molecule, second_molecule])
#         print('MCS atoms between 372 and 246 : ', mcs.numAtoms)
#         print('MCS SMARTS:', mcs.smartsString)

end = time.time()
print('Computed in: ', end-start)
# with Chem.SDMolSupplier('./DB00626.sdf') as suppl6:
#     ms6 = [x for x in suppl6 if x is not None]
#
# with Chem.SDMolSupplier('./DB00985.sdf') as suppl7:
#     ms7 = [x for x in suppl7 if x is not None]
#
# for molecule in ms:
#     for second_molecule in ms6:
#         print('Matches between molecules of 27 with 626: ', molecule.GetSubstructMatches(second_molecule, maxMatches=1000000))
#         print('Vice versa: ', second_molecule.GetSubstructMatches(molecule))
#         print(Chem.MolToSmiles(second_molecule))
#         print(Chem.MolToSmiles(molecule))
#         print('MCS atoms between 27 and 626 : ', rdFMCS.FindMCS([molecule, second_molecule]).numAtoms)
#
# for molecule in ms:
#     for second_molecule in ms7:
#         print('Matches between molecules of 27 with 985: ',
#               molecule.GetSubstructMatches(second_molecule, maxMatches=1000000))
#         print('Vice versa: ', second_molecule.GetSubstructMatches(molecule))
#         print(Chem.MolToSmiles(second_molecule))
#         print(Chem.MolToSmiles(molecule))
#         print('MCS atoms between 27 and 985 : ', rdFMCS.FindMCS([molecule, second_molecule]).numAtoms)
#

# for m in ms:
#     mol_id = m.GetProp('_Name')
#     smiles = Chem.MolToSmiles(m)
#     print(smiles)
#     # G.add_node(mol_id, smiles=smiles)
#
#     for atom in m.GetAtoms():
#         print(atom.GetSymbol(), atom)
#
# with Chem.SDMolSupplier('./DB00050.sdf') as suppl3:
#     ms3 = [x for x in suppl3 if x is not None]
#
# with Chem.SDMolSupplier('./DB00080.sdf') as suppl4:
#     ms4 = [x for x in suppl4 if x is not None]
#
# with Chem.SDMolSupplier('./DB00115.sdf') as suppl5:
#     ms5 = [x for x in suppl5 if x is not None]
#
# m = ms[0]
# m2 = ms2[0]
# m3 = ms3[0]
# m4 = ms4[0]
# m5 = ms5[0]
#
# print('Matches m with m2: ', m.GetSubstructMatches(m2))
# print('Matches m with m3: ', m.GetSubstructMatches(m3))
# print('Matches m with m4: ', m.GetSubstructMatches(m4))
# print('Matches m with m5: ', m.GetSubstructMatches(m5))
# print('Matches m with m: ', m.GetSubstructMatches(m))
#        # print(Chem.MolToSmiles(second_molecule))
        # print(Chem.MolToSmiles(molecule))
        # smarts = Chem.MolToSmarts(second_molecule)
        # print(smarts)
        # Chem.MolFromSmarts(smarts).ToSmiles()
        # print(Chem.MolToSmarts(molecule))
#