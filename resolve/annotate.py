import argparse
from cmath import exp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os, json
import seaborn as sns

# import scanpy as sc
import anndata

import napari
from skimage.io import imread
from skimage.transform import rescale
from dask import delayed
from anndata import AnnData
from napari.types import ImageData, PointsData, LabelsData, LayerDataTuple
from magicgui import magicgui
from magicgui import widgets
import pathlib

sns.set_style('ticks')
sns.set_context('paper')

from util import *

parser = argparse.ArgumentParser(description='Launch a GUI to manually annotate regions')
parser.add_argument('-d', '--datadir', required=True, help='The main data directory containing resolve data, ideally absolute path, e.g. [...]/data/')
parser.add_argument('-s', '--slide', default='C', help='The slide to plot, e.g. C or D')
parser.add_argument('-p', '--position', default='A1', help='The position on the slide to plot, A1 to D2')
parser.add_argument('-a', '--adata', help='The anndata file, for cell cluster visualization')
parser.add_argument('--experiment_id', default='P22344', help='ID of resolve experiment, for parsing file names')
parser.add_argument('--region_names', default='./region_names.csv', help='A file containing the hierarchical names of the regions to be annotated.')

def get_rgba(i):
    """
    Get rgba from a hex palette for cell clusters
    """
    palette=['#2f4f4f','#228b22','#00ff00','#000080','#1e90ff','#00ffff','#ff8c00','#deb887','#8b4513','#ff0000','#ff69b4','#800080']
    h = palette[int(i) % len(palette)].replace('#','')
    return [int(h[i:i+2], 16)/256 for i in (0, 2, 4)] + [0.8]

class AnnotationGUI(object):

    def __init__(self,
        adata, region_df,
        datadir='~/workspace/resolve/data',
        slide='M2', position='D1-1', experiment_id='P22344',
        ds=10):
        
        self.datadir = datadir
        self.slide_str = slide
        self.position_str = position
        self.ds_int = ds
        self.region_df = region_df

        self.points_range = (None, None)

        # self.dapi_array = self._read_dapi(
        #     os.path.join(datadir, f'{experiment_id}-slide{slide}/{experiment_id}-slide{slide}_{position}_DAPI.tiff'), ds=ds)
        self.gene_df = self._read_gene(
            os.path.join(datadir, f'segmentation_{experiment_id}-slide{slide}/{position}/segmentation.csv'))

        self.gene_list = adata.var.reset_index()['gene'].tolist()
        self.adata = adata
        
        region_anno_dir = os.path.join(datadir, 'region_annotation')
        if not os.path.isdir(region_anno_dir):
            os.makedirs(region_anno_dir)

    @staticmethod
    def _read_dapi(dapi_file, ds=10):
        """
        Read downsampled DAPI image to np array
        """
        lazy_imread = delayed(imread)
        reader = lazy_imread(dapi_file)  # doesn't actually read the file
        dapi = reader.compute()  # *now* it reads.
        
        return rescale(dapi, scale=1/ds, anti_aliasing=True)

    @staticmethod
    def _read_gene(gene_file):
        X = pd.read_table(gene_file, sep=',')
        X = X.query("~ is_noise")
        return X

    def _get_gene_exp(self, gene, linear_log='log'):
        is_in_slice = self.adata.obs.eval(f'slide == "{self.slide_str}" & slice == "{self.position_str}"')
        points_data = self.adata.obsm['spatial'][is_in_slice] / self.ds_int
        cell_size = np.sqrt(self.adata.obs[is_in_slice].area)
        cell_size = 25 * cell_size / np.max(cell_size)
        expression_level = self.adata.X[is_in_slice, self.adata.var.index == gene]
        if linear_log == 'linear':
            expression_level = (expression_level / np.max(expression_level))
        elif linear_log == 'log':
            expression_level = np.log(expression_level + 1)
            
        return dict(points_data=points_data, expression_level=expression_level, cell_size=cell_size)

    def run(self):
        """
        Launches GUI for annotation
        Args:
            adata - Anndata object
            datadir - main data dir containing resolve data, e.g. ./data/
            slide - str, 'C' or 'D'
            experiment_id - str, for parsing file names
        """

        @magicgui(auto_call=True,
                gene={"choices": self.gene_list, 'label': 'gene (per spot)'})
        def add_gene_layer(gene) -> LayerDataTuple:
            """
            Adds a slected gene to the plot
            Each point is a mRNA transcript
            The layer could be deleted from the GUI later
            """
            gene_palette = sns.color_palette('pastel')
            
            points_data = self.gene_df.query('gene == "%s"'%gene)[['x','y']].values / self.ds_int
            
            gene_color = gene_palette[np.random.choice(len(gene_palette))]
            gene_color = np.array(gene_color * points_data.shape[0]).reshape(-1,3)
            points_properties = {'name': gene, 
                                'face_color': gene_color, 
                                'size': 10}
            return points_data, points_properties, 'points'

        @magicgui(auto_call=True,
                gene={'choices': self.gene_list, 'label': 'gene (per cell)'},
                linear_log={'choices': ['linear', 'log']})
        def add_gene_expression_layer(gene, linear_log='log') -> LayerDataTuple:
            """
            Adds a slected gene to the plot
            Each point is a cell
            The layer could be deleted from the GUI later
            """
            gene_palette = sns.color_palette('pastel')
            
            gene_dict = self._get_gene_exp(gene, linear_log=linear_log)
            points_data = gene_dict['points_data']

            properties = {'expression_level': gene_dict['expression_level']}
            points_properties = {'name': gene,
                                'properties': properties,
                                'face_color': 'expression_level',
                                'face_colormap': 'inferno',
                                'size': gene_dict['cell_size'],
                                'edge_width': 0.0}
            return points_data, points_properties, 'points'

        @magicgui(call_button='Show coexpression',
                  geneA={'choices': self.gene_list, 'label': 'gene A (magenta)'},
                  geneB={'choices': self.gene_list, 'label': 'gene B (green)'},
                  cell_size_mode={'choices': ['fixed', 'actual'], 'label': 'cell size mode'})
        def add_gene_coexpression_layer(geneA, geneB, cell_size_mode='actual') -> LayerDataTuple:
            """
            Adds 2 genes and plot coexpression in 2 colors
            """
            def render_coexp_color(geneA_dict, geneB_dict):
                """
                Normalize **log transformed** expression level
                Returns:
                    face_color - (n_cell, 3) array, RGB
                """
                gene_coexp_colors = {'magenta': np.array([1,0,1]), 'green': np.array([0,1,0])}
                epsilon = 1e-3
                normexp = (geneA_dict['expression_level'] / np.maximum(np.percentile(geneA_dict['expression_level'], 99.9), epsilon),
                           geneB_dict['expression_level'] / np.maximum(np.percentile(geneB_dict['expression_level'], 99.9), epsilon))
                normexp = np.clip(normexp, 0, 1)
                face_color = (normexp[0].reshape(-1,1) * gene_coexp_colors['magenta'].reshape(1,3) +
                              normexp[1].reshape(-1,1) * gene_coexp_colors['green'].reshape(1,3))
                return face_color
                
            geneA_dict = self._get_gene_exp(geneA)
            geneB_dict = self._get_gene_exp(geneB)
            
            ## face_color - (n_cell, 3) array, RGB ##
            face_color = render_coexp_color(geneA_dict, geneB_dict)
            points_data = geneA_dict['points_data']
            if cell_size_mode == 'actual':
                cell_size = geneA_dict['cell_size']
            elif cell_size_mode == 'fixed':
                cell_size = 15
                
            points_properties = {'name': '%s & %s' % (geneA, geneB),
                    'face_color': face_color,
                    'size': cell_size,
                    'edge_width': 0.0}
            
            return points_data, points_properties, 'points'
            

        @magicgui(call_button='Save regions',
            filename={'widget_type': 'FileEdit'})
        def save_region(region_layer:LabelsData,
                        filename=os.path.join(self.datadir, 'region_annotation', f'slide{self.slide_str}_{self.position_str}')) -> None:
            np.save(filename, region_layer)


        @magicgui(call_button='Load regions',
            filename={'widget_type': 'FileEdit'})
        def load_region(filename=os.path.join(self.datadir, 'region_annotation', f'slide{self.slide_str}_{self.position_str}.npy')) -> LayerDataTuple:
            """
            Loads a npy file to the 'Regions' layer
            """
            try:
                region_data = np.load(filename)
            except:
                print('Not a .npy file')
            
            region_properties = {'name': 'Regions'}
            return region_data, region_properties
            

        mask = self.adata.obs['sample'] == 'slide%s_%s' % (self.slide_str, self.position_str)

        cell_coordinate = self.adata.obsm['spatial'].loc[mask,:] / self.ds_int
        x_shape, y_shape = int(np.ceil(np.max(cell_coordinate.x))), int(np.ceil(np.max(cell_coordinate.y)))
        cell_clusters = self.adata.obs.loc[mask,'clusters']
        cell_colors = [get_rgba(i) for i in cell_clusters]

        # viewer = napari.view_image(adjust_image(self.dapi_array))
        # cell_cluster_layer = viewer.add_points(cell_coordinate, symbol='square', face_color=cell_colors, name='Cell clusters')
        viewer = napari.view_image(0.1 * np.ones((x_shape, y_shape)), name='Background')
        cell_cluster_layer = viewer.add_points(cell_coordinate, symbol='square', size=10, face_color=cell_colors, name='Cell clusters')
        
        # labels_layer = viewer.add_labels(np.zeros_like(self.dapi_array, dtype=int), name='Regions')
        labels_layer = viewer.add_labels(np.zeros((x_shape, y_shape), dtype=int), name='Regions')

        labels_layer.features = self.region_df
        table_widget = widgets.Table(labels_layer.features)

        viewer.window.add_dock_widget(table_widget)
        viewer.window.add_dock_widget(add_gene_layer)
        viewer.window.add_dock_widget(add_gene_expression_layer)
        viewer.window.add_dock_widget(add_gene_coexpression_layer)
        viewer.window.add_dock_widget(save_region)
        viewer.window.add_dock_widget(load_region)

        # Shortcut to close the gui window
        @viewer.bind_key('q')
        def goodbye(viewer):
            print('Goodbye world!')
            viewer.close()
            
        # Hook up the points layer to the colorbar
        def print_layer_name(event):
            print(f"Active layer changed!")
            
        # viewer.layers.events.selection.events.active.connect(print_layer_name)

        # keep the dropdown menus in the gui in sync with the layer model
        viewer.layers.events.inserted.connect(add_gene_layer.reset_choices)
        viewer.layers.events.removed.connect(add_gene_layer.reset_choices)

        napari.run()


# @viewer.bind_key('b')
# def add_colorbar(viewer):
#     if viewer.points_range[0] is not None:
#         fig, ax = plt.subplots()
#         plt.colorbar(boundaries=self.points_range, ax=ax)
#         plt.show()

if __name__ == "__main__":

    args = parser.parse_args()

    adata_file = args.adata
    if adata_file is None:
        adata_file = os.path.join(args.datadir, 'adata', f'{args.experiment_id}_adata.h5ad')
    adata = anndata.read_h5ad(adata_file)
    region_df = pd.read_csv(args.region_names)

    gui = AnnotationGUI(
        adata = adata,
        datadir=args.datadir,
        slide=args.slide, position=args.position,
        region_df=region_df)

    gui.run()