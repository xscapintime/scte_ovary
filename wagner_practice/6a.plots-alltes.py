# %%
import logging, matplotlib, os, sys
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors

# %%
from glbase3 import *

# %%
plt.rcParams['figure.figsize']=(8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['font.size']=10

# %%
sc.settings.figdir = 'tes'

# %%
adata = sc.read('./learned.h5ad')
print(adata)
all_genes = adata.var['n_cells'].index # gene names are stored in the index

# %%
TEs = genelist(filename='../id_hg38/TE_genes_id.hg38.txt', format={'name': 0, 'force_tsv': True})

# %%
#%%time

for te in TEs:
    print(te['name'])
    if te['name'] in all_genes:
        sc.pl.umap(adata, color=[te['name'], te['name']], size=10, legend_loc='on data', show=False, save='TE-{0}.pdf'.format(te['name']), vmin=0, vmax=3)


