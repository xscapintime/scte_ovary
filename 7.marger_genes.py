import logging, matplotlib, os, sys
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
from rpy2.robjects.packages import importr
#from gprofiler import gprofiler
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 1
sc.set_figure_params(dpi=200, dpi_save=200)

sc.settings.figdir = 'markers'

adata = sc.read('learned.h5ad') # You can skip the script 3 if using te 2b.
#sc.pp.log1p(adata)

print(adata.var_names)

oh = open('gene_names.all.tsv', 'w')
for g in adata.var_names:
    oh.write('%s\n' % g)
oh.close()

marker_genes_dict = {
    'Spermatogonia': ['FGFR3', 'GAGE12J'],
    'Spermatocytes': ['DDX4', 'SYCP3',  'ZPBP'], # 'SPO11',
    'Spermatid': ['PRM3', 'LINC00467', 'TSACC'],
    'Sertoli': ['SOX9', 'MIR202HG',], # C4
    'Peritubular Myoid and Leydig': ['ACTA2', 'IGF1'],
    'Smooth Musc.': ['DCN'],
    'Mac': ['CD14', 'S100A4', 'FCER1G'],
    'Endothelial': ['VWF', 'GNG11', 'IFI27'],

    markers : ['UTF1', 'NODAL', 'GDF3', 'DPPA5', 'LEFTY1', 'LEFTY2', 'DPPA2',
    # HERVE-splice candidates:
    'SSUH2',
    'AC079466.1',
    'AC092666.1',
    'LINC01198',
    'SSUH2',
    'AP002495.1',
    ]

    }

sc.pl.stacked_violin(adata, marker_genes_dict, groupby='leiden_r0.10', rotation=90, dendrogram=True, show=False, save='markers.pdf')
sc.pl.dotplot(adata, marker_genes_dict, groupby='leiden_r0.10', dot_max=0.5, dendrogram=True, standard_scale='var', show=False, save='markers.pdf')
sc.pl.heatmap(adata, marker_genes_dict, groupby='leiden_r0.10', show=False, save='markers.pdf')

for k in marker_genes_dict:
    sc.pl.tsne(adata, color=marker_genes_dict[k], size=10, legend_loc='on data', vmax=3, show=False, save='markers-{0}.pdf'.format(k))
    sc.pl.umap(adata, color=marker_genes_dict[k], color_map='plasma', size=10, vmax=3, legend_loc='on data', show=False, save='markers-{0}.pdf'.format(k))

