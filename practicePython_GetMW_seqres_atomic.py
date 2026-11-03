from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import csv
import pandas as pd
import glob
import os



folder_path = r"D:/research v3.0/project associate/Complex Interface packing/1. Reading Manuscript 1 - protein packing/testing theory of tetramers and hexamers/structures_homoMW/extract"

df = []

# %% Reading the CIFs

for cif in glob.glob(os.path.join(folder_path, "*.cif")):   
    b = MMCIF2Dict(cif)
    file = os.path.basename(cif)
    print("reading "+file)
    # Collect all keys starting with "_entity." (excluding _entity_poly, _entity_src_gen, etc.)
    cols = [k for k in b if k.startswith("_entity.")]
    d2 = pd.DataFrame({k.replace("_entity.", ""): b[k] for k in cols})
    d2 = d2[d2['type'] =='polymer']
    mw_cif = 0
    comment = ""
    for i in range(len(d2)): 
        mw_cif = mw_cif + (float(d2.iat[i, 4])) * float((d2.iat[i, 5]))
        comment = comment + d2.iat[i, 4] + " x " + d2.iat[i, 5] + "; "
        
    
    
    
    # For reading seqres
    cols = [k for k in b if k.startswith("_pdbx_poly_seq_scheme.")]
    d2 = pd.DataFrame({k.replace("_pdbx_poly_seq_scheme.", ""): b[k] for k in cols})
    d2 = d2[d2['hetero'] =='n']

    d_seq = d2.drop(columns = ['ndb_seq_num','pdb_seq_num','auth_seq_num','pdb_mon_id','auth_mon_id','pdb_strand_id','pdb_ins_code'])
    d_atom = d2.drop(columns = ['seq_id','mon_id','ndb_seq_num','pdb_ins_code'])
    d_atom = d_atom[d_atom['pdb_mon_id'] != '?']

    ## now replacing res with weights


    cols = [k for k in b if k.startswith("_chem_comp.")]
    d2 = pd.DataFrame({k.replace("_chem_comp.", ""): b[k] for k in cols})
    d2 = d2[d2['type'] !='non-polymer']

    d_seq["weight"] = d_seq["mon_id"].map(d2.set_index("id")["formula_weight"])
    d_atom["weight"] = d_atom["pdb_mon_id"].map(d2.set_index("id")["formula_weight"])

    seq_len = len(d_seq)
    atom_len = len(d_atom)

    d_seq['weight'] = pd.to_numeric(d_seq['weight'], errors='coerce')
    mw_seq = d_seq['weight'].sum()

    d_atom['weight'] = pd.to_numeric(d_atom['weight'], errors='coerce')
    mw_atom = d_atom['weight'].sum()
    
    
    
    temp_data = { "file" : file, "mw_cif" : mw_cif,"mw_atom" : mw_atom, "mw_seq" : mw_seq, "seq_len" : seq_len, "atom_len" : atom_len, "comment" : comment }
    df.append(temp_data)
    
    
    
df = pd.DataFrame(df)    
df.to_csv("OwnCode_MWs_resLen.csv")
# %%test

#b = MMCIF2Dict("D:/research v3.0/project associate/Complex Interface packing/1. Reading Manuscript 1 - protein packing/testing theory of tetramers and hexamers/structures_homoMW/extract/1h5s-assembly1.cif")



# %%







#
#mw_cif = 0
#comment = ""
#for i in range(len(d2)): 
#    mw_cif = (d2.iat[i, 4]) * (d2.iat[i, 5])
#    comment = comment + d2.iat[i, 4] + " x " + d2.iat[i, 5] + "; "
    
#temp_data = { "file" : file, "mw_cif" : mw_cif, "comment" : comment }
#df.append(temp_data)
