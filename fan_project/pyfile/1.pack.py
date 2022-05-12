#!/usr/bin/env python
# coding: utf-8

# In[8]:


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


# In[9]:


samples = [
    sc_utils.sparsify("../fan2019_data/SRR7644615_merged.csv.gz",  obs_add={'sample_id': 'sample_1-1', 'ac': 'SRR7644615', 'source': 'stroma', 'patient': 'P9'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR7644616_merged.csv.gz",  obs_add={'sample_id': 'sample_1-2', 'ac': 'SRR7644616', 'source': 'stroma', 'patient': 'P7'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR7644618_merged.csv.gz",  obs_add={'sample_id': 'sample_1-3', 'ac': 'SRR7644618', 'source': 'stroma', 'patient': 'P7'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR7644619_merged.csv.gz",  obs_add={'sample_id': 'sample_1-4', 'ac': 'SRR7644619', 'source': 'stroma', 'patient': 'P9'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR7644621_merged.csv.gz",  obs_add={'sample_id': 'sample_1-5', 'ac': 'SRR7644621', 'source': '2-5mmFollicle', 'patient': 'P7'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR7644623_merged.csv.gz",  obs_add={'sample_id': 'sample_1-6', 'ac': 'SRR7644623', 'source': 'stroma', 'patient': 'P9'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR7644625_merged.csv.gz",  obs_add={'sample_id': 'sample_1-7', 'ac': 'SRR7644625', 'source': 'stroma', 'patient': 'P9'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR7644627_merged.csv.gz",  obs_add={'sample_id': 'sample_1-8', 'ac': 'SRR7644627', 'source': 'stroma', 'patient': 'P9'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR7644629_merged.csv.gz",  obs_add={'sample_id': 'sample_3-13', 'ac': 'SRR7644629', 'source': '1-2mmFollicle', 'patient': 'P3'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR7644631_merged.csv.gz",  obs_add={'sample_id': 'sample_3-14', 'ac': 'SRR7644631', 'source': '1-2mmFollicle', 'patient': 'P3'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR7644633_merged.csv.gz",  obs_add={'sample_id': 'sample_3-15', 'ac': 'SRR7644633', 'source': '1-2mmFollicle', 'patient': 'P3'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR7644635_merged.csv.gz",  obs_add={'sample_id': 'sample_3-16', 'ac': 'SRR7644635', 'source': '2-5mmFollicle', 'patient': 'P3'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR7644637_merged.csv.gz",  obs_add={'sample_id': 'sample_3-17', 'ac': 'SRR7644637', 'source': 'stroma', 'patient': 'P3'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR7644638_merged.csv.gz",  obs_add={'sample_id': 'sample_3-18', 'ac': 'SRR7644638', 'source': 'stroma', 'patient': 'P2'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR7644639_merged.csv.gz",  obs_add={'sample_id': 'sample_3-5', 'ac': 'SRR7644639', 'source': 'stroma', 'patient': 'P0'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR7644641_merged.csv.gz",  obs_add={'sample_id': 'sample_3-6', 'ac': 'SRR7644641', 'source': 'stroma', 'patient': 'P0'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR8428465_merged.csv.gz",  obs_add={'sample_id': 'sample1_B1_i12A', 'ac': 'SRR8428465', 'source': 'stroma', 'patient': 'P7'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR8428466_merged.csv.gz",  obs_add={'sample_id': 'sample10_B2_i10D', 'ac': 'SRR8428466', 'source': '1-2mmFollicle', 'patient': 'P7'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR8428467_merged.csv.gz",  obs_add={'sample_id': 'sample11_B2_i10F', 'ac': 'SRR8428467', 'source': 'stroma', 'patient': 'P9'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR8428468_merged.csv.gz",  obs_add={'sample_id': 'sample12_B2_i10G', 'ac': 'SRR8428468', 'source': 'stroma', 'patient': 'P0'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR8428469_merged.csv.gz",  obs_add={'sample_id': 'sample13_B1_i12C', 'ac': 'SRR8428469', 'source': 'stroma', 'patient': 'P0'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR8428470_merged.csv.gz",  obs_add={'sample_id': 'sample145_B1_i12D', 'ac': 'SRR8428470', 'source': '2-5mmFollicle', 'patient': 'P7'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR8428471_merged.csv.gz",  obs_add={'sample_id': 'sample2_B1_i12B', 'ac': 'SRR8428471', 'source': 'stroma', 'patient': 'P7'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR8428472_merged.csv.gz",  obs_add={'sample_id': 'sample3_B1_i12G', 'ac': 'SRR8428472', 'source': 'stroma', 'patient': 'P7'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR8428473_merged.csv.gz",  obs_add={'sample_id': 'sample4_B1_i12F', 'ac': 'SRR8428473', 'source': 'stroma', 'patient': 'P7'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR8428474_merged.csv.gz",  obs_add={'sample_id': 'sample5_B2_i10E', 'ac': 'SRR8428474', 'source': 'stroma', 'patient': 'P7'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR8428475_merged.csv.gz",  obs_add={'sample_id': 'sample6a_B1_i12H', 'ac': 'SRR8428475', 'source': 'stroma', 'patient': 'P7'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR8428476_merged.csv.gz",  obs_add={'sample_id': 'sample7_B2_i10C', 'ac': 'SRR8428476', 'source': '1-2mmFollicle', 'patient': 'P7'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR8428477_merged.csv.gz",  obs_add={'sample_id': 'sample8a_B2_i10A', 'ac': 'SRR8428477', 'source': '2-5mmFollicle', 'patient': 'P7'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR8428478_merged.csv.gz",  obs_add={'sample_id': 'sample8b_B2_i10B', 'ac': 'SRR8428478', 'source': '2-5mmFollicle', 'patient': 'P7'}, csv=True),
    sc_utils.sparsify("../fan2019_data/SRR8428479_merged.csv.gz",  obs_add={'sample_id': 'sampleC1_B1_i12E', 'ac': 'SRR8428479', 'source': '2-5mmFollicle', 'patient': 'P7'}, csv=True)
]
print('Loaded Samples...')


# In[10]:


# Quick pre-filtering, these should be low, otherwise it can mess up downstream analysis, but also can get rid of trivial uninteresting things
[sc.pp.filter_cells(sam, min_genes=100) for sam in samples]
[sc.pp.filter_cells(sam, max_counts=200000) for sam in samples]
[sc.pp.filter_cells(sam, min_counts=100) for sam in samples]

# Do not filter gene here; concatenate joins on the union, so if a gene fails in a single sample, it will also be deleted from all other samples;


# In[11]:


print('Concatenating')
adata = samples[0].concatenate(samples[1:])


# In[12]:


del samples


# In[13]:


adata.X = adata.X.astype('float32')


# In[14]:


print(adata)


# In[15]:


print('Total number of cells: {:d}'.format(adata.n_obs))
print('Total number of genes: {:d}'.format(adata.n_vars))


# In[16]:


adata.write('./raw_data.h5ad')


# In[17]:


oh = open('gene_names.all.tsv', 'w')
for g in adata.var_names:
    oh.write('%sn' % g)
oh.close()
