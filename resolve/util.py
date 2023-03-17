import numpy as np
import pandas as pd
# import scanpy as sc
import matplotlib.pyplot as plt
import os, json
import seaborn as sns
# import anndata
# from anndata import AnnData
import scipy.stats as stats
from skimage.io import imread
from tqdm import tqdm
from matplotlib.backends.backend_pdf import PdfPages
from anndata import AnnData
tqdm.pandas()
sns.set_style('ticks')
sns.set_context('paper')

def save_fig(filename, fig=None):

    figdir, _ = os.path.split(filename)
    if not os.path.isdir(figdir):
        os.makedirs(figdir)

    if fig is None:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
    else:
        fig.savefig(filename, dpi=300, bbox_inches='tight')


def save_multi_image(filename):
    pp = PdfPages(filename)
    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()

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

def adjust_image(image, down_sample=1):
    return np.rot90(np.fliplr(image[::down_sample, ::down_sample]), k=1)


def plot_slide(samples, adata_all, slide2plot='C', figname=None):
    slide_samples = samples.query('slide == "%s"'%slide2plot)#['sliceID']

    fig, ax = plt.subplots(4, 2, figsize=(12,20))
    
    for i,letter in enumerate(['A', 'B', 'C', 'D']):
        for j,number in enumerate([1, 2]):
            position = '%s%d'%(letter, number)

            if position in slide_samples['slice'].values:
                
                _ = sc.pl.embedding(adata_all[adata_all.obs['sample'] == 'slide%s_%s' % (slide2plot, position),:],
                                     basis="spatial", color="clusters", ax=ax[i,j], show=False)
                ax[i,j].set_title(position)
    plt.suptitle('slide%s'%slide2plot)
    
    if figname is None:
        figname = './fig/slide%s_spatial.pdf'%slide2plot
    save_fig(figname, fig)


def count_gene_slide(samples, gene, adata_all, slide='C', use_raw=True, percentile=None, cutoff=None, figname=None):
    slide_samples = samples.query('slide == "%s"'%slide)#['sliceID']

    if percentile is not None:
        if use_raw:
            cutoff = np.percentile(adata_all[adata_all.obs['slide'] == slide, adata_all.var.index == gene].layers['raw'], percentile)
        else:
            cutoff = np.percentile(adata_all[adata_all.obs['slide'] == slide, adata_all.var.index == gene].X, percentile)
            
        print('cutoff:', cutoff)
    elif cutoff is None:
        raise ValueError('Either percentile or absolute cutoff should be given')
    
    fig, ax = plt.subplots(4, 2, figsize=(12,20))
    
    for i,letter in enumerate(['A', 'B', 'C', 'D']):
        for j,number in enumerate([1, 2]):
            position = '%s%d'%(letter, number)

            if position in slide_samples['slice'].values:
                adata = adata_all[adata_all.obs['sample'] == 'slide%s_%s' % (slide, position),:]
                
                gene_adata = adata[:, adata.var.index == gene]
                
                if use_raw:
                    gene_counts = gene_adata.layers['raw'].flatten()
                else:
                    gene_counts = gene_adata.X.flatten()
                
                n_positive = np.sum(gene_counts > cutoff)
                n_total = len(gene_counts)
                
                sns.ecdfplot(gene_counts, color='k', ax=ax[i,j])
                ax[i,j].plot([0,np.max(gene_counts)], np.array([.01,.01])*percentile, 'r--')

                ax[i,j].set_title('%s: %d/%d, %.2f%%' % (position, n_positive, n_total, 100*n_positive/n_total))
                
    plt.suptitle('slide%s, cutoff=%.2f' % (slide, cutoff))
    
    if figname is None:
        figname = './fig/slide%s_%s_cell_count.pdf' % (slide, gene)
    save_fig(figname, fig)


def read_cell_data(cell, segmentation_dir_prefix='segmentation'):
    """
    Returns cropped DAPI, raw and molecule arrays from cell name,
    eg. '3031-slideC_C1'
    """
    def get_ind(coor, max_coor):
        crop_ind = min(max(0, coor), max_coor)
        return int(crop_ind)
    
    m = 100

    cell_id = int(cell.split('-')[0])
    slide, tile = cell.split('-')[1].split('_')

    segmentation_dir = './data/%s_32801-%s/%s'%(segmentation_dir_prefix, slide, tile)
    spot_file = os.path.join(segmentation_dir, 'segmentation.csv')
    cell_stats_file = os.path.join(segmentation_dir, 'segmentation_cell_stats.csv')
    raw_data_dir = './data/32801-%s'%slide
    dapi_file = os.path.join(raw_data_dir, '32801-%s_%s_DAPI.tiff'%(slide, tile))
    raw_file = os.path.join(raw_data_dir, '32801-%s_%s_raw.tiff'%(slide, tile))

    cell_stats = pd.read_csv(cell_stats_file, index_col='cell')
    cell_stat = cell_stats[cell_stats.index == cell_id]
    spot = pd.read_csv(spot_file, sep=',').query('cell == %d'%cell_id)
    
    # print(xmin, xmax, ymin, ymax)
    dapi = imread(dapi_file)
    imx, imy = dapi.shape
    xmin, xmax, ymin, ymax = get_ind(spot.y.min()-m, imx), get_ind(spot.y.max()+m, imx), get_ind(spot.x.min()-m, imy), get_ind(spot.x.max()+m, imy)
    dapi = dapi[xmin:xmax, ymin:ymax]
    raw = imread(raw_file)[xmin:xmax, ymin:ymax]
    
    spots = np.zeros_like(dapi)
    for i,s in spot.iterrows():
        try:
            spots[int(s.y - xmin), int(s.x - ymin)] += 1
        except:
            print('crop indices: ', xmin, xmax, ymin, ymax)
            print(imx, imy)
        
    return cell_stat, spot, spots, dapi, raw, (xmin, xmax, ymin, ymax)


def plot_cell(cell, segmentation_dir_prefix='segmentation'):
    cell_stat, spot, spots, dapi, raw, crop_ind = read_cell_data(cell, segmentation_dir_prefix)
    fig, ax = plt.subplots(1, 3, figsize=(12,4), sharey=True, sharex=True)
    # sns.scatterplot(y=spot.x-crop_ind[2], x=spot.y-crop_ind[0], hue=spot.gene, palette='tab20', ax=ax[0])
    # if len(np.unique(spot.gene)) > 10:
    #     ax[0].get_legend().remove()
    ax[0].imshow(np.rot90(spots), vmax=0.2, cmap='gray')
    ax[0].set_title('%s, area = %d'%(cell, cell_stat.area.values[0]))
    ax[1].imshow(np.rot90(dapi), cmap='gray')
    ax[2].imshow(np.rot90(raw), cmap='gray')

"""
Functions to help load region segmentation .npy files from the gui
back into anndata
"""

def get_region_annotation(spatial, annotation, region_names):
    """
    Apply region annotation output from the GUI to adata.obs
    via coordinates in adata.obsm['spatial']
    """
    X, Y = annotation.shape
    ds = 10
    idx = int(np.clip(spatial.x/ds, 0, X-1))
    idy = int(np.clip(spatial.y/ds, 0, Y-1))
    region_id = annotation[idx, idy]
    return region_names.loc[region_id,:]

def add_region_annotation_to_adata(adata, annotation_file, region_names_file, sample_name=None):
    """
    Args:
        adata - Anndata object
        annotation_file - e.g. './resolve_analysis/resolve/slideC_A1.npy'
        regions_names_file - e.g. './resolve_analysis/resolve/region_names.csv'
        sample - str, name of the sample to annotate to. Matches adata.obs.sample
            !!! If not provided, assume all cells are on the same slide (dangerous)
    Returns:
        adata - Anndata object with updated adata.obs field
    """
    # adata_all = anndata.read_h5ad('./data/32801_resolve_adata.h5ad')
    annotation = np.load(annotation_file)
    region_names = pd.read_csv(region_names_file)
    
    if sample_name is None:
        adata_slide = adata.copy()
    else:
        adata_slide = adata[adata.obs['sample'] == sample_name]
        
    cell_annotation = adata_slide.obsm['spatial'].progress_apply(lambda row: get_region_annotation(row, annotation, region_names), axis=1)
    adata_slide.obs = adata_slide.obs.join(cell_annotation, how='left')

    return adata_slide
    
"""
*may be useless
anndata helper functions
"""
def get_gene_idx(adata, gene):
    if not gene in adata.var.index:
        print("gene %s not in adata" % gene)
        return np.nan
        
    idx = np.argmax(adata.var.index == gene)
    return idx
    
def get_gex_vec(adata, gene):
    idx = get_gene_idx(adata, gene)
    return adata.X[:,idx]