import pandas as pd
import general_loader
import yaml
import sqlite3

from Bio import SeqIO
from Bio.KEGG import REST
# from Bio.KEGG.KGML import KGML_parser
# from Bio.Graphics.KGML_vis import KGMLCanvas
from time import sleep
from random import randint
import random
import io
import os

import datetime
import signal

import pickle

import debugpy as debug
debug.listen(("0.0.0.0", 5678))
print("waiting attach!")
debug.wait_for_client()

def update(conn, task):
    sql = ''' UPDATE compound_pert
              SET query = 1 
              WHERE rowid = ?'''
    cur = conn.cursor()
    r = cur.executemany(sql, task)
    conn.commit()
    return cur.rowcount

# def to_df(result):
#     return pd.read_table(io.StringIO(result), header=None)

def to_df(result):
    try:
        r = pd.read_table(io.StringIO(result), sep='delimiter', header=None ,engine='python')
    except Exception as e:
        r = pd.DataFrame({})
    return r 

#长时间运行信号
def signal_handler(signum, frame):
    raise Exception("KEGG API Timed out!")

#网络查询KEGG
#每个cpd都会得到一些drug数据
def query_network():
    conn = sqlite3.connect('l1000.db')
    #GEO92742 50K 
    #GSE70139 2K 
    compounds = pd.read_sql('''
                                select rowid,pert_id,pert_iname,pert_type,query 
                                from compound_pert
                                where query = '0'
                            ''',conn)

    df1 = pd.DataFrame()
    temporary_compounds = []
    print("Query KEGG")

    from tqdm import tqdm
    for i in tqdm(range(len(compounds))):
        query_result = '\n'
        compound = compounds.iloc[i]
        # query_result = REST.kegg_find("drug", compound['pert_iname']).read()
        temporary_compounds.append(compound['rowid'].astype(str))

        try:
            signal.signal(signal.SIGALRM, signal_handler)
            signal.alarm(10)   # 10 秒
            query_result = REST.kegg_find("drug", compound['pert_iname']).read()
            
        except Exception as e:
            now = datetime.datetime.now()
            with open('error.log', 'a') as file:
                file.write('{t};{id};{c};{e}\n'.format(t=now,id=compound['pert_id'],c=compound['pert_iname'],e=str(e)))
                
        #临时
        if query_result != '\n':
                df2 = to_df(query_result)
                df2['pert_id'] = compound['pert_id']
                df1 = pd.concat([df1, df2], ignore_index = True).drop_duplicates()
        
        if i % 10 == 0 and i > 0 and len(df1) > 0:
            df1.to_csv("query_lincs_compound_kegg.csv", mode='a')
            df1 = pd.DataFrame()
            #q = '\''+'\',\''.join(temporary_compounds)+'\''
            q = (((val,) for val in temporary_compounds))
            r = update(conn,q)
            temporary_compounds = []
            sleep(random.uniform(0.1, 1.9))
    df1.to_csv("query_lincs_compound_kegg.csv", mode='a')
    conn.close()

# query_network()
#======================处理数据========================

def count():
    conn = sqlite3.connect('l1000.db')
    sigs = pd.read_sql_query('''
                              select cs.cell_id,cd.*
                              from compound_sig cs
                              left join compound_data cd on cs.sig_id = cd.sig_id
                              where cs.pert_id = 'BRD-K57080016'
                              ''',conn)
    sigs_sum_data = (sigs.iloc[:,2:]).sum(axis=1, numeric_only=True)
    sig_cells = sigs.iloc[:,:1]
    result = pd.concat([sig_cells,sigs_sum_data ], axis=1)
    result.columns = ['cell_id','sum_value']
    print("end count")
#count()


def handle_kegg_query_result(kegg_result,result_dic):
    G = kegg_result[kegg_result['pert_id']!= 'pert_id']
    n = ['index','drug_id','name','pert_id']
    G.columns = n
    G = G.groupby(['pert_id'],as_index=False).sum()
    for i in range(len(G)):
        pert_id = G.iloc[i]['pert_id']
        if pert_id not in result_dic:
            result_dic[pert_id] = []
        else:
            print(pert_id)
        drug_ids = G.iloc[i]['drug_id'].split(':')
        for drugs_per_cpd in drug_ids:
            val = drugs_per_cpd.replace('dr','')
            if val != '':
                result_dic[pert_id].append(val)
                
    return result_dic

#网络查询结束后，得到query_lincs_compound_kegg.csv
#处理该csv.转为字典类型:cpd:{drugs},得到keqq_query_result.pkl
def process_kegg_query_format_into_dict():
    kegg_GSE = pd.read_csv('query_lincs_compound_kegg.csv')
    result = {}
    result = handle_kegg_query_result(kegg_GSE,result)
    
    with open('keqq_query_result.pkl', 'wb') as handle:
        pickle.dump(result, handle)
        
# process_kegg_query_format_into_dict()

# 利用 keqq_query_result.pkl 中处理好的cpd-drugs
# 查询 drugs 对应的
def create_raw_data():
    start_time = datetime.datetime.now()
    with open('keqq_query_result.pkl', 'rb') as f: 
        query_result = pickle.load(f)
    dict = {}
    try:
        with open('ground_truth.pkl', 'rb') as loaded: 
            loaded_file = pickle.load(loaded)
    except Exception as e:
        loaded_file = {}
        
    diff=list(set(query_result.keys()).difference(loaded_file.keys()))
    #len(diff),len(loaded_file.keys()),len(query_result.keys())
    #(1681, 1987, 3668)
    #(1552, 2116, 3668)
    #(1433, 2235, 3668)
    #(1201, 2467, 3668)
    for key in diff:
        drugs = query_result[key] # key = cmap
        key_dict = {}
        print(key,datetime.datetime.now())
        for drug in drugs:
            a = ''
            try:
                signal.signal(signal.SIGALRM, signal_handler)
                signal.alarm(20)   # 20 秒
                a = REST.kegg_get(drug).read()
            except Exception as e:
                now = datetime.datetime.now()
                with open('ground_truth_error.log', 'a') as file:
                    file.write('{t};{id};{c};{e}\n'.format(t=now,id=key,c=drug,e=str(e)))
            #a = REST.kegg_get(drug).read()
            b = to_df(a)
            pathway = []
            start = 0
            for i in range(len(b)):
                row = b.iloc[i].values[0]
                d = []
                if 'PATHWAY' in row:
                    start = i
                    d = row.replace('PATHWAY',' ').split('  ')
                    while('' in d):
                        d.remove('')
                    pathway.append(d[1])
                elif start > 0 and row.index(' ') == 0 if (' ' in row) else False:
                    d = row.split('  ')
                    while('' in d):
                        d.remove('')
                    pathway.append(d[1])
                elif start > 0 and row.index(' ') > 0 if (' ' in row) else False:
                    break; 
            if len(pathway) > 0: 
                key_dict[drug] = pathway
        try:
            with open('ground_truth.pkl','rb') as saved:
                dict = pickle.load(saved)
                dict[key] = key_dict
        except:
            dict[key] = key_dict
        with open('raw_data.pkl', 'wb') as handle:
            pickle.dump(dict, handle)
#create_raw_data()

def create_ground_truth():
    import pandas as pd
    with open('raw_data.pkl', 'rb') as f: 
        raw_data = pickle.load(f)
    ground_truth_col = []
    
    for cpd in raw_data.keys():
        cpd_raw = raw_data[cpd]
        cpd_pathways = []
        for drug in cpd_raw.keys():
            pathways = cpd_raw[drug]
            cpd_pathways = cpd_pathways + pathways
        ground_truth_col = list(set(ground_truth_col + cpd_pathways))
    ground_truth_col = ['pert_name'] + ground_truth_col    
    df = pd.DataFrame(columns = ground_truth_col)
    for cpd in raw_data.keys():
        cpd_raw = raw_data[cpd]
        cpd_data_dict = {'pert_name':cpd}
        for drug in cpd_raw.keys():
            pathways = cpd_raw[drug]
            for pathway in pathways:
                cpd_data_dict[pathway] = 1
        cpd_data = pd.DataFrame([cpd_data_dict])
        df = pd.concat([df,cpd_data],ignore_index=True)
    print(df)
    df_data = df.fillna(0)
    df_data.to_csv('ground_truth.csv')
#create_ground_truth()

print('done')