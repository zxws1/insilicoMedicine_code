import torch
from torch.utils.data import DataLoader
from L1000DataSet import L1000DataSet
#from sklearn.metrics.pairwise import cosine_similarity
from torch.nn.functional import cosine_similarity
from model import JointEncoder
import debugpy as debug
# debug.listen(("0.0.0.0", 5678))
# print("waiting attach!")
# debug.wait_for_client()
def train():
	PATH = './a.pth'
	model = JointEncoder()
	data = L1000DataSet()
	batch_size = 16
	prepared_data = DataLoader(dataset=data, batch_size=batch_size, shuffle=True, drop_last=True)
	print(len(data),batch_size)
	learning_rate = 1e-4
	optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=1e-8)
	cos_similar = torch.nn.CosineSimilarity(dim=0, eps=1e-8)
	for i, batch_value in enumerate(prepared_data):
		ecfps_p = batch_value[0]
		ecfps_n = batch_value[1]
		genes = batch_value[2]
		#print(ecfps.shape,genes.shape)
		p,n,a = model(ecfps_p.float(),ecfps_n.float(),genes.float())
		loss = 1/torch.sum(torch.minimum((cos_similar(p,a) - cos_similar(n,a) - torch.abs(p-n)),torch.tensor(0)))
		print(loss.item())
		optimizer.zero_grad()
		loss.backward()
		optimizer.step()
	torch.save(model.state_dict(), PATH)

	print("end")
	data.close_conn()
train()
