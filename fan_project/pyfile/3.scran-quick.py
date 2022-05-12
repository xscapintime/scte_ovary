#!/usr/bin/env python
# coding: utf-8

# In[1]:


import logging, os, sys
import rpy2
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import anndata2ri
import anndata
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
#import rpy2.rinterface_lib.callbacks
import anndata2ri
from matplotlib import rcParams
from matplotlib import colors


# In[2]:


plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=300)
sc.settings.figdir = 'markers'
sc.settings.autoshow = False


# In[3]:


# Ignore R warning messages
#Note: this can be commented out to get more verbose R output
#rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)
# Automatically convert rpy2 outputs to pandas dataframes
pandas2ri.activate()
anndata2ri.activate()
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
#sc.set_figure_params(dpi=200, dpi_save=300)
sc.logging.print_versions()


# In[4]:


adata = sc.read('filtered.h5ad')
print(adata)
print('Number of cells: {:d}'.format(adata.n_obs))


# In[5]:


adata_pp = adata.copy()  # not used later?
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=20, svd_solver='arpack')
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added='groups', resolution=0.3)


# In[6]:


input_groups = adata_pp.obs['groups']
data_mat = adata.X.T
data_mat = sp.sparse.csc_matrix(data_mat)



# In[10]:


# Skip this if it's the second time for time-saving
#from scipy.io import mmwrite
#mmwrite('sparse_data_mat.mtx', data_mat)


# In[11]:


#get_ipython().magic('load_ext rpy2.ipython')


# In[12]:


#get_ipython().run_cell_magic('R', '', '\nlibrary(Matrix)\nlibrary(scran)')


# In[13]:


#get_ipython().run_cell_magic('R', '', "\ndata_mat_r <- readMM('sparse_data_mat.mtx')\nstr(data_mat_r)")


# In[ ]:


#get_ipython().run_cell_magic('time', '', '\n%%R  -i input_groups -o size_factors\n\nsize_factors = computeSumFactors(data_mat_r, clusters=input_groups, min.mean=0.1)')


# In[ ]:

#np.save("size_factors.npy", size_factors)


# In[18]:


del adata_pp


# In[19]:

size_factors = np.load("size_factors.npy")
adata.obs['size_factors'] = size_factors


# In[20]:


adata.layers['counts'] = adata.X.copy()
adata.X /= adata.obs['size_factors'].values[:,None]
sc.pp.log1p(adata)
adata.X = np.clip(adata.X, 0, 100000000) # get rid of <0
adata.X = sp.sparse.csr_matrix(adata.X)
adata.raw = adata # You only need to do this if you do batch correction


# In[ ]:


# combat batch correction could go here:
sc.pp.combat(adata, key='sample_id') # by sample itself?


# In[ ]:


# resparsify:
adata.X = np.clip(adata.X, 0, 1e9) # get rid of <0
adata.X = sp.sparse.csr_matrix(adata.X)

adata.write('./normed.h5ad')
