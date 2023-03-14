import numpy as np
import pandas as pd
# import scanpy as sc
import matplotlib.pyplot as plt
import os, json
import seaborn as sns
# import anndata
from anndata import AnnData
import scipy.stats as stats
# from skimage.io import imread
from tqdm import tqdm


def read_adata(segmentation_dir):
    """
    Returns adata object from segmentation output dir of
    a single tile, e.g. A1
    """
    count_file = os.path.join(segmentation_dir, 'segmentation_counts.tsv')
    cell_stats_file = os.path.join(segmentation_dir, 'segmentation_cell_stats.csv')
    
    anno = pd.read_csv(cell_stats_file, index_col='cell')
    anno.index = anno.index.astype(str)
    X = pd.read_csv(count_file, sep='\t', index_col='gene').T
    adata = AnnData(X=X,
                    obs=anno.iloc[:,2:], layers={'raw': X})
    adata.obsm['spatial'] = anno.loc[:,['x','y']]
    return adata