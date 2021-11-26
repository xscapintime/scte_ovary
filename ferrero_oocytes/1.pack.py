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
    sc_utils.sparsify("ferrero2019/ss.SRR8446786.1.csv.gz", obs_add={'sample': 'DS10_S1', 'age': '31YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446787.1.csv.gz", obs_add={'sample': 'DS11_S2', 'age': '31YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446788.1.csv.gz", obs_add={'sample': 'DS12_S3', 'age': '31YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446789.1.csv.gz", obs_add={'sample': 'DS13_S4', 'age': '31YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446790.1.csv.gz", obs_add={'sample': 'DS14_S5', 'age': '24YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446791.1.csv.gz", obs_add={'sample': 'DS15_S6', 'age': '24YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446792.1.csv.gz", obs_add={'sample': 'DS16_S7', 'age': '24YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446793.1.csv.gz", obs_add={'sample': 'DS17_S8', 'age': '24YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446796.1.csv.gz", obs_add={'sample': 'DS4_S4', 'age': '22YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446797.1.csv.gz", obs_add={'sample': 'DS5_S5', 'age': '22YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446798.1.csv.gz", obs_add={'sample': 'DS6_S6', 'age': '22YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446799.1.csv.gz", obs_add={'sample': 'DS7_S7', 'age': '27YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446800.1.csv.gz", obs_add={'sample': 'DS1_S1', 'age': '18YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446801.1.csv.gz", obs_add={'sample': 'DS2_S17', 'age': '18YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446802.1.csv.gz", obs_add={'sample': 'DS2_S2', 'age': '18YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446803.1.csv.gz", obs_add={'sample': 'DS3_S3', 'age': '18YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446804.1.csv.gz", obs_add={'sample': 'DS8_S8', 'age': '27YO'}, csv=True),
    sc_utils.sparsify("ferrero2019/ss.SRR8446805.1.csv.gz", obs_add={'sample': 'DS9_S9', 'age': '27YO'}, csv=True),
    ]
print('Loaded Samples...')


# %%
## SMART-seq, for cell for one sample, fake filtering mito-genes, for calculatring merics
# for sam in samples:
#     sam.var['mt'] = sam.var_names.str.startswith('mt-')
#     sam.var['mt']


# [sc.pp.calculate_qc_metrics(sam, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True) for sam in samples]


# %%
# Quick pre-filtering, these should be low, otherwise it can mess up downstream analysis, but also can get rid of trivial uninteresting things
[sc.pp.filter_cells(sam, min_genes=1000) for sam in samples]
# [sc.pp.filter_cells(sam, max_counts=200000) for sam in samples]
[sc.pp.filter_cells(sam, min_counts=5000) for sam in samples]


# %%
# plot for each sample
for sam in samples:
    sc.pl.violin(
        sam, 
        [
         'n_genes', 
         'n_counts'
         ],
        multi_panel=True,
        save='qc1-pre-samples_{:s}'.format(samples[0].obs.index[0]) + '.pdf'
    )


# Do not filter gene here; concatenate joins on the union, so if a gene fails in a single sample, it will also be deleted from all other samples;
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
