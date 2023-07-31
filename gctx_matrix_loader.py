import os
import dask.array as da
import dask.dataframe as dd
import h5py
import pandas as pd

# HDF5 nodes
ROW_META_GROUP_NODE = "/0/META/ROW"
COL_META_GROUP_NODE = "/0/META/COL"
RID_NODE = "/0/META/ROW/id"
CID_NODE = "/0/META/COL/id"
DATA_NODE = "/0/DATA/0/matrix"



class GctxMatrixLoader:
    def __init__(self, path, sample_index_name, gene_index_name):
        self.path = path
        # self.sample_index_name = sample_index_name
        # self.gene_index_name = gene_index_name
        # self.n_rows, self.n_cols = self.get_shape()

    @classmethod
    def from_node(cls, node, data_dir):
        path = os.path.join(data_dir, node.pop("file"))
        return cls(path, **node)


    def read(self, sample_ids=None, max_genes=None):
        """Returns an expression dataframe subsetted by held out samples (row_index)"""
        with h5py.File(self.path, "r") as gctx_file:
            sample_ids = pd.DataFrame(gctx_file[CID_NODE]).astype(str)
            gene_ids = pd.DataFrame(gctx_file[RID_NODE]).astype(str)
            data = pd.DataFrame(gctx_file[DATA_NODE])
            data.index = sample_ids[0]
            data.columns = gene_ids[0]
        return data,sample_ids,gene_ids