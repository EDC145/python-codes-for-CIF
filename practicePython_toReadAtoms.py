# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 10:29:19 2026

@author: Eshaan
"""

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import csv
import pandas as pd
import glob
import os



folder_path = r"D:/research v3.0/project associate/Complex Interface packing/1. Reading Manuscript 1 - protein packing/testing theory of tetramers and hexamers/structures_homoMW/extract"

df = []


# %%
for cif in glob.glob(os.path.join(folder_path, "*.cif")):   
    b = MMCIF2Dict(cif)
    file = os.path.basename(cif)
    print("reading "+file)
    
    # For reading atomic coordinates
    cols = [k for k in b if k.startswith("_atom_site.")]
    d2 = pd.DataFrame({k.replace("_atom_site.", ""): b[k] for k in cols})
    d2 = d2[d2['label_atom_id'] == 'CA']
    d2 = d2[d2['Cartn_x'] != '?']

    true_atom = d2

    cols = [k for k in b if k.startswith("_chem_comp.")]
    d2 = pd.DataFrame({k.replace("_chem_comp.", ""): b[k] for k in cols})
    d2 = d2[d2['type'] !='non-polymer']
    true_atom["weight"] = true_atom["label_comp_id"].map(d2.set_index("id")["formula_weight"])

    true_atom['weight'] = pd.to_numeric(true_atom['weight'], errors='coerce')
    mw_atom = true_atom['weight'].sum()


    temp_data = { "file" : file, "mw_atom" : mw_atom }
    df.append(temp_data)
    
    
    
df = pd.DataFrame(df)    
df.to_csv("D:/research v3.0/project associate/Complex Interface packing/5. correcting homomer sets/my own code/OwnCode_MWs_fromCoordinates_2.csv")




# %% test

b =  MMCIF2Dict("D:/research v3.0/project associate/Complex Interface packing/1. Reading Manuscript 1 - protein packing/testing theory of tetramers and hexamers/structures_homoMW/extract/6qv4-assembly1.cif")
