#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:


plt.rcParams['figure.figsize'] = (10,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 6
sc.settings.autoshow = False


# In[3]:


adata = sc.read('raw_data.h5ad')


# In[4]:


#mito_genes = adata.var_names.str.startswith('mt-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
#adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
#adata.obs['n_counts'] = adata.X.sum(axis=1).A1


# In[6]:


## violin plot
sc.pl.violin(adata, ['n_genes', 'n_counts'], groupby='source', size=0, log=False, cut=0, show=False, save='qc1-pre-norm-source_v2.pdf')
sc.pl.violin(adata, ['n_genes', 'n_counts'], groupby='patient', size=0, log=False, cut=0, show=False, save='qc1-pre-norm-patient_v2.pdf')


# In[7]:


# Base filtering for QC failures:
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=4000)
sc.pp.filter_cells(adata, min_counts=500)
sc.pp.filter_cells(adata, max_counts=20000)
sc.pp.filter_genes(adata, min_cells=3) # Only filter genes here; ## default value, also used in fan et al's paper
#adata = adata[adata.obs['percent_mito'] < 0.2, :]


# In[8]:


# Filter out the mt genes to stop them etting into the most variable
#mito_genes = adata.var_names.str.startswith('mt-')
#mask = np.isin(adata.var_names, mito_genes, invert=True, assume_unique=True)
#adata = adata[:, mask] # also slices .var[]
# The above causes trouble for some reason in scran


# In[10]:


## violin plot after QC filtering
sc.pl.violin(adata, ['n_genes','n_counts'], groupby='source', size=0, log=False, cut=0, show=False, save='qc1-source_v2.pdf')
sc.pl.violin(adata, ['n_genes','n_counts'], groupby='patient', size=0, log=False, cut=0, show=False, save='qc1-patient_v2.pdf')


# In[11]:


print('Total number of cells: {:d}'.format(adata.n_obs))
print('Total number of genes: {:d}'.format(adata.n_vars))


# In[12]:


adata.write('./filtered.h5ad')


# In[13]:


oh = open('gene_names.filtered.tsv', 'w')
for g in adata.var_names:
    oh.write('%s\n' % g)
oh.close()
