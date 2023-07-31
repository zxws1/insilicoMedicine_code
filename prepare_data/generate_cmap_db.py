import pandas as pd
import yaml
import debugpy as debug
import csv_metadata_loader
import sqlite3
import gctx_matrix_loader
import pickle

# debug.listen(("0.0.0.0", 5678))
# print("waiting attach!")
# debug.wait_for_client()

GP_yaml_file = 'GP.yaml'

def load_mate(yaml_file='',mate='pert_metadata'):
    with open(yaml_file) as f:
        d = yaml.safe_load(f)
    data_dir = d["data_dir"]
    loader = csv_metadata_loader.CsvMetaDataLoader.from_node(d[mate], data_dir)
    df = loader.read()
    return df

def load_CMAP_data_to_sqlite():
    conn = sqlite3.connect('cmap.db')
    sig_meta = 'signature_metadata'
    df_GP_sig = load_mate(yaml_file=GP_yaml_file,mate=sig_meta)
    df_GP_sig.to_sql("GP_sig",conn,if_exists='replace',index=False)
    
    gene_meta = 'gene_metadata'
    df_GP_gene = load_mate(yaml_file=GP_yaml_file,mate=gene_meta)
    df_GP_gene.to_sql('GP_gene',conn,if_exists='replace',index=False)
    
    pert_meta = 'pert_metadata'
    df_GP_compound = load_mate(yaml_file=GP_yaml_file,mate=pert_meta)
    df_GP_compound.to_sql('GP_compound',conn,if_exists='replace',index=False)
    conn.commit()
    conn.close()
    
    with open(GP_yaml_file) as f:
        d = yaml.safe_load(f)
    data_dir = "/data/datacenter/H3C_GPU/projects/zhouxing/GP"
    #xpr
    xpr = {'sample_index_name': 'sig_id', 'gene_index_name': 'gene_id', 'file':'xpr.gctx'}
    data_loader_xpr = gctx_matrix_loader.GctxMatrixLoader.from_node(xpr, data_dir)
    #cp
    cp = {'sample_index_name': 'sig_id', 'gene_index_name': 'gene_id', 'file':'cp.gctx'}
    data_loader_cp = gctx_matrix_loader.GctxMatrixLoader.from_node(cp, data_dir)
    #oe
    oe = {'sample_index_name': 'sig_id', 'gene_index_name': 'gene_id', 'file':'oe.gctx'}
    data_loader_oe = gctx_matrix_loader.GctxMatrixLoader.from_node(oe, data_dir)
    #sh
    sh = {'sample_index_name': 'sig_id', 'gene_index_name': 'gene_id', 'file':'sh.gctx'}
    data_loader_sh = gctx_matrix_loader.GctxMatrixLoader.from_node(sh, data_dir)
    
    GP_data_xpr,GP_sample_ids_xpr,GP_gene_ids_xpr =  data_loader_xpr.read()
    GP_data_cp, GP_sample_ids_cp, GP_gene_ids_cp = data_loader_cp.read()
    GP_data_oe,GP_sample_ids_oe,GP_gene_ids_oe =  data_loader_oe.read()
    GP_data_sh, GP_sample_ids_sh, GP_gene_ids_sh = data_loader_sh.read()
    
    conn = sqlite3.connect('l1000.db')
    df_GP_key = pd.read_sql("select * from gp_gene where feature_space = 'landmark'",conn)
    key_genes = list(map(str, df_GP_key[['gene_id']].values.flatten().tolist()))
    key_names = list(map(str, df_GP_key[['gene_symbol']].values.flatten().tolist()))
    conn.close()
    GP_data_xpr = GP_data_xpr.loc[:, key_genes]
    GP_data_cp = GP_data_cp.loc[:, key_genes]
    GP_data_oe = GP_data_oe.loc[:,key_genes]
    GP_data_sh = GP_data_sh.loc[:,key_genes]
    
    frames = [GP_data_xpr, GP_data_cp, GP_data_oe, GP_data_sh]
    result = pd.concat(frames)

    result.columns = key_names
    result.reset_index(inplace=True)
    
    conn = sqlite3.connect('cmap.db')
    result.to_sql("GP_data",conn,if_exists='replace',index=False)
    
    conn.commit()
    conn.close()
    # 执行 ALTER TABLE GP_data RENAME COLUMN "0" TO "sig_id"
    
#load_CMAP_data_to_sqlite()

#获取与gene相关的signature，用于predication
def Get_CMAP_signature_by_genes():
    import datetime
    conn = sqlite3.connect('cmap.db')
    #CMAP是cpd
    sql = '''select gd.* 
            from GP_data gd 
            left join GP_sig gs on gd.sig_id = gs.sig_id 
            where gs.cmap_name in ( select gene_symbol from GP_gene )
            '''
    start = datetime.datetime.now()
    gp_data = pd.read_sql(sql,conn)
    end = datetime.datetime.now()
    print((end-start).total_seconds())
    with open('GP_embedding.pkl', 'wb') as handle:
        pickle.dump(gp_data, handle)
Get_CMAP_signature_by_genes()
print('done')
