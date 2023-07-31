import pandas as pd
import general_loader

import debugpy as debug
debug.listen(("0.0.0.0", 5678))
print("waiting attach!")
debug.wait_for_client()

#====================GP========================
GP_yaml_file = 'GP.yaml'
GP_data_loader, GP_sample_meta_loader, GP_gene_meta_loader = general_loader.set_loader(GP_yaml_file)

df_GP_samples = GP_sample_meta_loader.read()
#df_GP_genes = GP_gene_meta_loader.read()
GP_data,GP_sample_ids,GP_gene_ids = GP_data_loader.read()

a = df_GP_samples.join(GP_data,how='inner')
a.to_pickle("./GP.pkl")

#====================compound========================
compound_yaml_file = 'compound.yaml'
compound_data_loader, compound_sample_meta_loader, compound_gene_meta_loader = general_loader.set_loader(compound_yaml_file)
df_compound_samples = compound_sample_meta_loader.read()
#df_compound_genes = compound_gene_meta_loader.read()
compound_data,compound_sample_ids,compound_gene_ids = compound_data_loader.read()

b = df_compound_samples.join(compound_data,how='inner')
b.to_pickle("./compound.pkl")

# #====================================================
# import csv_metadata_loader
# import yaml
# with open(compound_yaml_file) as f:
#     d = yaml.safe_load(f)
# data_dir = d["data_dir"]
# GSM_compound_metadata_loader = csv_metadata_loader.CsvMetaDataLoader.from_node(d["pert_metadata"], data_dir)
# GSM_df_compound = GSM_compound_metadata_loader.read()

# with open(GP_yaml_file) as f:
#     d = yaml.safe_load(f)
# data_dir = d["data_dir"]
# GP_compound_metadata_loader = csv_metadata_loader.CsvMetaDataLoader.from_node(d["compound_metadata"], data_dir)
# GP_df_compound = GP_compound_metadata_loader.read()
print('done')