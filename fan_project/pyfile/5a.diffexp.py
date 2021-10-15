#!/usr/bin/env python
# coding: utf-8

# In[1]:


import logging, matplotlib, os, sys, glob
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import pandas as pd
from glbase3 import genelist, glload


# In[2]:


plt.rcParams['figure.figsize']=(8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['font.size']=10



# In[4]:


sc.settings.figdir = 'diffexp'

# In[]:
def mkdir(path):
 
    folder = os.path.exists(path)
 
    if not folder:                   #判断是否存在文件夹如果不存在则创建为文件夹
        os.makedirs(path)            #makedirs 创建文件时如果路径不存在会创建这个路径
        print("---  new folder...  ---")
        print("---  OK  ---")
 
    else:
        print("---  There is this folder!  ---")
        
gls = "/mnt/e/projects/scte-gpc/scte_ovary/fan_project/pyfile/gls"
mkdir(gls)

# In[5]:


[os.remove(f) for f in glob.glob('{0}/*.pdf'.format(sc.settings.figdir))]
[os.remove(f) for f in glob.glob('gls/*.glb'.format(sc.settings.figdir))]
[os.remove(f) for f in glob.glob('gls/*.tsv'.format(sc.settings.figdir))]


# In[6]:


de_leiden = 'leiden_r0.10'


# In[7]:


adata = sc.read('./learned.h5ad')
sc.tl.rank_genes_groups(adata, de_leiden, method='t-test_overestim_var', n_genes=3000)
#sc.tl.filter_rank_genes_groups(adata, min_fold_change=1)
adata.write('./de.h5ad')