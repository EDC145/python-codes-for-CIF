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



folder_path = r"D:/research v3.0/project associate/Complex Interface packing/7. Bulk Packing attempt - 1/selected_CIFs"
output_path = r"D:/research v3.0/project associate/Complex Interface packing/7. Bulk Packing attempt - 1/Only_Atomic_Coords/"
df = []


# %%
for cif in glob.glob(os.path.join(folder_path, "*.cif")):   
    b = MMCIF2Dict(cif)
    file = os.path.basename(cif)
    print("reading "+file)
    
    # For reading atomic coordinates
    cols = [k for k in b if k.startswith("_atom_site.")]
    d2 = pd.DataFrame({k.replace("_atom_site.", ""): b[k] for k in cols})
 #   d2 = d2[d2['label_atom_id'] == 'CA']
    d2 = d2[d2['Cartn_x'] != '?']
    d2 = d2[d2['group_PDB'] == 'ATOM']
    d2 = d2[d2['label_alt_id'].isin(['.', 'A'])] ## Need to check if others also
    d2 = d2[d2['label_comp_id'].isin(['ALA', 'ASP', 'THR', 'ILE', 'VAL', 'GLU', 'LEU', 'TYR', 'PRO', 'ASN', 'GLY','CYS',   'SER', 'HIS', 'LYS', 'ARG', 'TRP', 'MET', 'GLN', 'PHE'])]
    d2 = d2[d2['label_atom_id'].isin(['N',  'CA', 'SG' , 'C',   'O',  'CB',  'CG', 'OD1', 'OD2', 'OG1', 'CG2', 'CG1',     'CD1',  'CD', 'OE1', 'OE2', 'CD2', 'CE1', 'CE2',  'CZ',  'OH', 'ND2',  'OG','ND1', 'NE2',  'CE',  'NZ',  'NE', 'NH1', 'NH2', 'NE1', 'CE3', 'CZ2', 'CZ3','CH2',  'SD', 'OXT','AD1', 'AD2', 'AE1', 'AE2'])]    
    d2['occupancy'] = pd.to_numeric(d2['occupancy'], errors='coerce')
    d2 = d2[d2['occupancy'] > 0]
    models = d2['pdbx_PDB_model_num'].unique()      # in order of appearance
    if len(models) > 1:
        print(f"WARNING {file}: {len(models)} models found, keeping model {models[0]}")
        # to skip such files instead of keeping model 1, use: continue
    d2 = d2[d2['pdbx_PDB_model_num'] == models[0]]
    
    
    d3 =d2.drop(columns=['group_PDB','type_symbol','label_atom_id','label_alt_id',    'label_comp_id',                    'pdbx_PDB_ins_code','occupancy',                        'B_iso_or_equiv', 'Cartn_x_esd', 'Cartn_y_esd', 'Cartn_z_esd',                       'occupancy_esd', 'B_iso_or_equiv_esd', 'pdbx_formal_charge',                   'auth_atom_id',   'auth_seq_id','auth_comp_id','auth_asym_id'              ,       'pdbx_PDB_model_num','label_asym_id', 'label_entity_id', 'label_seq_id'])
    d3['id'] = range(1, len(d3) + 1)
    df = pd.DataFrame(d3)
    out_file = os.path.join(output_path, file + "_coord.txt")
    print(f"Writing file to {out_file}")
    d3.to_csv(out_file, sep=" ", header=False, index=False)
    
    
    
    
#df = pd.DataFrame(df)    
#df.to_csv("D:/research v3.0/project associate/Complex Interface packing/5. correcting homomer sets/my own code/OwnCode_MWs_fromCoordinates_2.csv")

#df.to_csv(file + "_coord_fromCIF", sep="\t", header=False, index=False)


# %% test

b =  MMCIF2Dict("D:/research v3.0/project associate/Complex Interface packing/7. Bulk Packing attempt - 1/selected_CIFs/6ahg-assembly1.cif")
