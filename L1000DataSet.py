from torch.utils.data import Dataset
from torch.utils.data import DataLoader
import torch
import pandas as pd
import sqlite3
from rdkit import Chem
from rdkit.Chem import AllChem,DataStructs
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import random

# import debugpy as debug
# debug.listen(("0.0.0.0", 5678))
# print("waiting attach!")
# debug.wait_for_client()

def smiles2fp(smilesstr,sig_id=0,pert_row_id=0):
	mol = Chem.MolFromSmiles(smilesstr)
	fp_obj = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048,
												   useChirality=True)
	arr = np.zeros((1,))
	DataStructs.ConvertToNumpyArray(fp_obj, arr)

	return arr

class L1000DataSet(Dataset):
	def __init__(self):
		self.data_path = 'l1000.db'
		self.genes_list = torch.range(1, 113871,dtype=torch.long)# list(range(118050))
		self.compound_list = torch.range(1, 1797,dtype=torch.long)# list(range(2170))
		self.conn = sqlite3.connect(self.data_path)
	def __len__(self):
		genes_ids = (pd.read_sql('''
                            select cd.rowid
                            from compound_data cd 
							left join compound_sig cs on cd.sig_id = cs.sig_id
							left join compound_pert cp on cs.pert_id = cp.pert_id
							where cp.canonical_smiles is not null and cp.canonical_smiles != 'restricted'
                           ''',self.conn))
		self.genes_list = genes_ids[['rowid']].values.flatten()
		compound_ids = (pd.read_sql("select rowid from compound_pert",self.conn))
		self.compound_list = compound_ids[['rowid']].values.flatten()
		return len(self.genes_list)

	def __getitem__(self, idx):
		index = idx+1
  
		sig_rowid = self.genes_list[index-1] # rowid [1,max index + 1]
		matrix_data = (pd.read_sql("select * from compound_data where rowid = {sig}".format(sig=sig_rowid),self.conn))
		genes = matrix_data.iloc[:,1:] #anchor
		sig_id = matrix_data[['sig_id']].values.flatten()[0]
		sig = pd.read_sql("select * from compound_sig where sig_id = '{sig}'".format(sig=sig_id),self.conn)
		pert_id = sig[['pert_id']].values.flatten()[0]
		positive_compound = pd.read_sql("select rowid,* from compound_pert where pert_id = '{pert}'".format(pert=pert_id),self.conn)
		positive_smiles = positive_compound[['canonical_smiles']].values.flatten()[0]
  
		selected_compound_id = positive_compound[['rowid']].values.flatten()[0]
		other_compound_index = random.randint(1,len(self.compound_list))
		while(other_compound_index == selected_compound_id):
			other_compound_index = random.randint(1,len(self.compound_list))
		other_compound_rowid = self.compound_list[other_compound_index-1]# rowid [1,max index + 1]
		negative_compound = pd.read_sql("select rowid,* from compound_pert where rowid = '{pert}'".format(pert=other_compound_rowid),self.conn)
		negative_smiles = negative_compound[['canonical_smiles']].values.flatten()[0]
  
		positive_smiles_ecfp = torch.from_numpy(smiles2fp(positive_smiles,sig_id=sig_id))
		negative_smiles_ecfp = torch.from_numpy(smiles2fp(negative_smiles,pert_row_id={selected_compound_id:other_compound_index}))
		return positive_smiles_ecfp,negative_smiles_ecfp,torch.tensor(genes.values.astype(float)).squeeze(0)
		
	def close_conn(self):
		self.conn.close()
		pass
def test():
	data = L1000DataSet()
	batch_size = 2
	prepared_data = DataLoader(dataset=data, batch_size=batch_size, shuffle=True, drop_last=True)
	print(len(data),batch_size)
	for i, batch_value in enumerate(prepared_data):
		ecfps = batch_value[0]
		genes = batch_value[1]
		print(ecfps.shape,genes.shape)
		#cosine_similarity()
		if i == 4:
			data.close_conn()
			break
	print("end")