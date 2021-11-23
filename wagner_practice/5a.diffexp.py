# %%
import logging, matplotlib, os, sys, glob
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import pandas as pd
from glbase3 import genelist

# %%
plt.rcParams['figure.figsize']=(8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['font.size']=10

# %%
from glbase3 import genelist, glload

# %%
sc.settings.figdir = 'diffexp'

# %%
[os.remove(f) for f in glob.glob('{0}/*.pdf'.format(sc.settings.figdir))]
[os.remove(f) for f in glob.glob('gls/*.glb'.format(sc.settings.figdir))]
[os.remove(f) for f in glob.glob('gls/*.tsv'.format(sc.settings.figdir))]

# %%
de_leiden = 'leiden_r0.10'

# %%
adata = sc.read('./learned.h5ad')
sc.tl.rank_genes_groups(adata, de_leiden, method='t-test_overestim_var', n_genes=3000)
#sc.tl.filter_rank_genes_groups(adata, min_fold_change=1)
adata.write('./de.h5ad')


