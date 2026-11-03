from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import csv
import pandas as pd
import glob
import os



a = pd.read_csv("D:/research v3.0/project associate/Complex Interface packing/5. correcting homomer sets/my own code/MW_output_chgpt_mod.csv")



folder_path = r"D:/research v3.0/project associate/Complex Interface packing/1. Reading Manuscript 1 - protein packing/testing theory of tetramers and hexamers/structures_homoMW/extract"

newWeight = []

# %%


for cif in glob.glob(os.path.join(folder_path, "*.cif")):    
    b = MMCIF2Dict(cif)
    file = os.path.basename(cif)
    print("reading "+file)
    max_mol_weight = max([float(w) for w in b['_entity.formula_weight'] if w not in ('?', '.')])
    print("weight is " + str(max_mol_weight))
    temp_data = { "file" : file, "max_mol_weight" : max_mol_weight }
    newWeight.append(temp_data)
    print("Added to dataframe")


# %%
newWeight = pd.DataFrame(newWeight)

#a = a[1:len(a)]

a.set_index('File', inplace = True)
newWeight.set_index('file', inplace = True)
c = a.join(newWeight, lsuffix='_left', rsuffix='_right', how='outer').reset_index()

c.to_csv("Comparing MW CIF.csv")
