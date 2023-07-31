import yaml
import csv_metadata_loader
import gctx_matrix_loader

def set_loader(yaml_file):
    with open(yaml_file) as f:
        d = yaml.safe_load(f)
    data_dir = d["data_dir"]
    data_loader = gctx_matrix_loader.GctxMatrixLoader.from_node(d["data"], data_dir)
    sample_meta_loader = csv_metadata_loader.CsvMetaDataLoader.from_node(d["sample_metadata"], data_dir)
    gene_meta_loader = csv_metadata_loader.CsvMetaDataLoader.from_node(d["gene_metadata"], data_dir)
    return data_loader, sample_meta_loader, gene_meta_loader

def set_lincs_loader(yaml_file):
    with open(yaml_file) as f:
        d = yaml.safe_load(f)
    data_dir = d["data_dir"]
    data_loader = gctx_matrix_loader.GctxMatrixLoader.from_node(d["data"], data_dir)
    sig_meta_loader = csv_metadata_loader.CsvMetaDataLoader.from_node(d["signature_metadata"], data_dir)
    pert_meta_loader = csv_metadata_loader.CsvMetaDataLoader.from_node(d["pert_metadata"], data_dir)
    return data_loader, sig_meta_loader, pert_meta_loader