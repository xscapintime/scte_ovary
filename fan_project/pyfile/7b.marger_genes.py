#!/usr/bin/env python
# coding: utf-8

# In[2]:


# %load 7b.marger_genes.py
import logging, matplotlib, os, sys
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
from glbase3 import glload
#from gprofiler import gprofiler
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 1
sc.set_figure_params(dpi=200, dpi_save=200)

sc.settings.figdir = 'herve_spliced'


# In[3]:


adata = sc.read('learned.h5ad') # You can skip the script 3 if using te 2b.

print(adata.var_names)

oh = open('gene_names.all.tsv', 'w')
for g in adata.var_names:
    oh.write('%s\n' % g)
oh.close()


# In[4]:


genes_to_do =[
    # HERVE-splice candidates:
    'SSUH2',
    'AC079466.1',
    'AC092666.1',
    'LINC01198',
    'SSUH2',
    'AP002495.1',
    ]


# In[5]:


for gene_name in genes_to_do:
    if gene_name not in adata.var_names:
        print('{0} gene name not found in lookup!'.format(gene_name))
        continue

    try:
        sc.pl.umap(adata, color=gene_name,
        size=10,
        legend_loc='on data',
        vmax=3,
        show=False,
        save='markers-{0}.pdf'.format(gene_name)
        )

    except KeyError: # this specific transcript_id is missing;
        print('{0} {1} not found'.format(transcript['transcript_id'], transcript['name']))


# In[ ]:




