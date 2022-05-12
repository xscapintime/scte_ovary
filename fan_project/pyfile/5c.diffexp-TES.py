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


plt.rcParams['figure.figsize']=(4,4)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['font.size']=10


# In[4]:


sc.settings.figdir = 'diffexp-tes'


# In[5]:


[os.remove(f) for f in glob.glob('{0}/*.pdf'.format(sc.settings.figdir))]
[os.remove(f) for f in glob.glob('gls/*.glb'.format(sc.settings.figdir))]
[os.remove(f) for f in glob.glob('gls/*.tsv'.format(sc.settings.figdir))]


# In[6]:


de_leiden = 'leiden_r0.10'


# In[7]:


adata = sc.read('./de.h5ad')
TEs = set(genelist(filename='../../id_hg38/TE_genes_id.hg38.txt', format={'name': 0, 'force_tsv': True})['name'])


# In[8]:


sc.pl.rank_genes_groups(adata, n_genes=25, sharey=True, show=False, save='genes-top25.pdf')
sc.pl.rank_genes_groups(adata, key='rank_genes_groups', show=False, save='genes.pdf')
sc.pl.rank_genes_groups_dotplot(adata, key='rank_genes_groups', show=False, save='genes-top25.pdf')


# In[9]:


#print(pd.DataFrame(adata.uns['rank_genes_groups']))
print(pd.DataFrame(adata.uns['rank_genes_groups']['names']))


# In[10]:


print()
topall = pd.DataFrame(adata.uns['rank_genes_groups']['names']) # get all;
fcs = pd.DataFrame(adata.uns['rank_genes_groups']['logfoldchanges'])
padj = pd.DataFrame(adata.uns['rank_genes_groups']['pvals_adj'])


# In[11]:


get_ipython().run_cell_magic('time', '', "\n# Matrix of DE genes:\n\ngroups = list(topall.columns.values)\n\nnewcols = {}\n\nfor group in groups:\n    newcols[group] = []\n\n    t = zip([i[group] for i in adata.uns['rank_genes_groups']['names']], [i[group] for i in adata.uns['rank_genes_groups']['logfoldchanges']], [i[group] for i in adata.uns['rank_genes_groups']['pvals_adj']])\n\n    print('Group: {0}'.format(group))\n    print(t)\n\n    for item in t:\n        if item[0] not in TEs:\n            continue\n\n        if item[1] < 2: # fold change\n            continue\n        if item[2] > 0.01: # just in case\n            continue\n\n        newcols[group].append({'name': item[0], 'log2FC': item[1], 'q': item[2]})")


# In[12]:


get_ipython().run_cell_magic('time', '', "\n# join all and draw a dotplot:\nfor group in newcols:\n    print(newcols[group])\n    if newcols[group]:\n        gl = genelist()\n        gl.load_list(newcols[group])\n        gl.saveTSV('gls/de_genes-grp{0}.tsv'.format(group))\n        gl.save('gls/de_genes-grp{0}.glb'.format(group))\n\n        genes = [i['name'] for i in newcols[group]]\n        sc.pl.dotplot(adata, genes, groupby=de_leiden, dot_max=0.7, dendrogram=True, standard_scale='var', show=False, save='de-grp{0}.pdf'.format(group))\n        sc.pl.matrixplot(adata, genes, groupby=de_leiden, dendrogram=True, standard_scale='var', show=False, save='de-grp{0}.pdf'.format(group))")


# In[13]:


get_ipython().run_cell_magic('time', '', "\nfor grp in newcols:\n    if not newcols[grp]:\n        continue\n    for k in newcols[grp]:\n        title = k['name']\n        sc.pl.umap(adata, color=k['name'], size=20, legend_loc='on data',\n            title=title,\n            vmin=0, vmax=3,\n            show=False, save='-markers-grp{0}-{1}.pdf'.format(grp, k['name']))\n        #sc.pl.violin(adata, [k], groupby='disease', size=0, log=False, cut=0, show=False, save='markers-{0}-disease.pdf'.format(k))\n        #sc.pl.violin(adata, [k], groupby='cell_type', size=0, log=False, cut=0, show=False, save='markers-{0}-cell_type.pdf'.format(k))")
