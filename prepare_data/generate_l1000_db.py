import pandas as pd
import general_loader
import yaml
import sqlite3
import csv_metadata_loader

# import debugpy as debug
# debug.listen(("0.0.0.0", 5678))
# print("waiting attach!")
# debug.wait_for_client()

def load_mate(yaml_file='',mate='compound_metadata'):
    with open(yaml_file) as f:
        d = yaml.safe_load(f)
    data_dir = d["data_dir"]
    loader = csv_metadata_loader.CsvMetaDataLoader.from_node(d[mate], data_dir)
    df = loader.read()
    return df

def Load_LINCS_to_sqlite():
    #关键基因
    conn = sqlite3.connect('l1000.db')
    gp_yaml_file = 'GP.yaml'
    gene_mate = 'gene_metadata'
    df_GP_gene = load_mate(yaml_file=gp_yaml_file,mate=gene_mate)
    r = df_GP_gene.to_sql("gp_gene",conn,if_exists='replace',index=False)
    conn.commit()
    conn.close()

    conn = sqlite3.connect('l1000.db')
    df_GP_key = pd.read_sql("select * from gp_gene where feature_space = 'landmark'",conn)
    key_genes = list(map(str, df_GP_key[['gene_id']].values.flatten().tolist()))
    
    yaml_file_70138 = 'compound_70138.yaml'
    data_loader_70138, sig_meta_loader_70138, pert_meta_loader_70138 = general_loader.set_lincs_loader(yaml_file_70138)
    df_sig_70138 = sig_meta_loader_70138.read()
    df_pert_70138 = pert_meta_loader_70138.read()
    data_70138,sample_ids_70138,gene_ids_70138 = data_loader_70138.read()
    data_70138 = data_70138.loc[:, key_genes]
    
    yaml_file_92742 = 'compound_92742.yaml'
    data_loader_92742, sig_meta_loader_92742, pert_meta_loader_92742 = general_loader.set_lincs_loader(yaml_file_92742)
    df_sig_92742 = sig_meta_loader_92742.read()
    df_pert_92742 = pert_meta_loader_92742.read()
    data_92742,sample_ids_92742,gene_ids_92742 = data_loader_92742.read()
    data_92742 = data_92742.loc[:, key_genes]
    
    sig = [df_sig_70138, df_sig_92742]
    sig_df = pd.concat(sig)
    sig_df.reset_index(inplace=True)
    
    pert = [df_pert_70138, df_pert_92742]
    pert_df = pd.concat(pert)
    pert_df.reset_index(inplace=True)
    
    data = [data_70138, data_92742]
    data_df = pd.concat(data)
    data_df.reset_index(inplace=True)
    
    sig_df.to_sql("compound_sig",conn,if_exists='replace',index=False)
    pert_df.to_sql("compound_pert",conn,if_exists='replace',index=False)
    data_df.to_sql("compound_data",conn,if_exists='replace',index=False)
    
    conn.commit()
    conn.close()
    # 执行 ALTER TABLE compound_pert ADD query CHAR(25) DEFAULT '0'

Load_LINCS_to_sqlite()
print("done")
