# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
"""

Pack the scRNA-seq data using scanpy, prep for scran normalisation

"""

import logging, matplotlib, os, sys
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from anndata import AnnData
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
from rpy2.robjects.packages import importr
plt.rcParams['figure.figsize'] = (8, 8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 10
sc.settings.autoshow = False

import sc_utils


# %%
samples = [
    sc_utils.sparsify("wagner2020/ss.C-sec-1.csv.gz",  obs_add={'replicate': "C-sec#1", 'type': 'C-sec'}, csv=True),
    sc_utils.sparsify("wagner2020/ss.C-sec-2.csv.gz",  obs_add={'replicate': "C-sec#2", 'type': 'C-sec'}, csv=True),
    sc_utils.sparsify("wagner2020/ss.GRP-1.csv.gz", obs_add={'replicate': "GRP#1", 'type': 'GRP'}, csv=True),
    sc_utils.sparsify("wagner2020/ss.GRP-2.csv.gz", obs_add={'replicate': "GRP#2", 'type': 'GRP'}, csv=True),
    ]
print('Loaded Samples...')


# %%
# Quick pre-filtering, these should be low, otherwise it can mess up downstream analysis, but also can get rid of trivial uninteresting things

[sc.pp.filter_cells(sam, min_genes=1000) for sam in samples]
[sc.pp.filter_cells(sam, max_counts=200000) for sam in samples]
[sc.pp.filter_cells(sam, min_counts=5000) for sam in samples]

# Do not filter gene here; concatenate joins on the union, so if a gene fails in a single sample, it will also be deleted from all other samples;


# %%
print('Concatenating')
adata = samples[0].concatenate(samples[1:])


# %%
del samples


# %%
adata.X = adata.X.astype('float32')


# %%
print(adata)


# %%
print('Total number of cells: {:d}'.format(adata.n_obs))
print('Total number of genes: {:d}'.format(adata.n_vars))


# %%
adata.write('./raw_data.h5ad')


# %%
oh = open('gene_names.all.tsv', 'w')
for g in adata.var_names:
    oh.write('%s\n' % g)
oh.close()


