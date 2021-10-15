#!/usr/bin/env python
# coding: utf-8

# In[1]:


import logging, matplotlib, os, sys
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors


# In[2]:


from glbase3 import *


# In[3]:


plt.rcParams['figure.figsize']=(8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['font.size']=10


# In[4]:


sc.settings.figdir = 'tes'


# In[5]:


adata = sc.read('./learned.h5ad')
print(adata)
all_genes = adata.var['n_cells'].index # gene names are stored in the index


# In[6]:


TEs = genelist(filename='./id_hg38/TE_genes_id.hg38.txt', format={'name': 0, 'force_tsv': True})


# In[7]:


get_ipython().run_cell_magic('time', '', "\nfor te in TEs:\n    print(te['name'])\n    if te['name'] in all_genes:\n        sc.pl.umap(adata, color=[te['name'], te['name']], size=10, legend_loc='on data', show=False, save='TE-{0}.pdf'.format(te['name']), vmin=0, vmax=3)")


# In[ ]:





# In[ ]:





# In[ ]:




