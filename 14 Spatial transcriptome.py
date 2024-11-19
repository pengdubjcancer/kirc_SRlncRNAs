import os
import anndata as ad
import squidpy as sq
import pandas as pd
import scanpy as sc
import numpy as np
from scipy.sparse import csr_matrix
import matplotlib.pylab as plt
import seaborn as sns 
from umap.umap_ import UMAP
from ALLCools.clustering import tsne, significant_pc_test, log_scale
from ALLCools.plot import 
import warnings
warnings.filterwarnings("ignore")
sc.settings.set_figure_params(dpi=800)
os.getcwd()
os.chdir("/home/ccrccsp")

adata = sc.read_h5ad("942d1f67-b0be-4177-84f9-90689c534b42.h5ad")
adata = sc.read_h5ad("aae49bb1-9dc3-4dd7-88a5-7c6cf7f17948.h5ad")
adata = sc.read_h5ad("dd27114b-6e1b-4fd7-b5f7-489f5d79dcf8.h5ad") 

adata = sc.read_h5ad("6800STDY12499406.h5ad") 
adata
a = adata.var
b = adata.obs
print(b.columns.values)
c = adata.uns['spatial']
d = adata.uns["spatial"]["spaceranger130_count_44317_6800STDY12499406_GRCh38-2020-A"]
sc.tl.pca(adata) 
sc.pp.neighbors(adata)
sc.tl.umap(adata)
# sc.tl.leiden(adata) 
del c["is_single"]
del adata.uns['spatial']
adata.uns['spatial'] = c 
sc.pl.spatial(adata, img_key="hires", color=['ENSG00000256128'], size=1.5,color_map='Reds',save = "6800STDY12499406LINC00944.pdf") 
sc.pl.spatial(adata, img_key="hires", color=['ENSG00000229847'], size=1.5,color_map='Reds',save = "6800STDY12499406EMX2OS.pdf") 
sc.pl.spatial(adata, img_key="hires", color="cell_type", size=1.5,save = "6800STDY12499406celltype.pdf")
sc.pl.spatial(adata, img_key="hires", size=0,save = "6800STDY12499406he.pdf") 
