#!/usr/bin/env python
# coding: utf-8

# In[2]:


# %load 7.marger_genes.py
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


# In[3]:


adata = sc.read('learned.h5ad') # You can skip the script 3 if using te 2b.
#sc.pp.log1p(adata)

print(adata.var_names)

oh = open('gene_names.all.tsv', 'w')
for g in adata.var_names:
    oh.write('%s\n' % g)
oh.close()


# In[5]:


marker_genes_dict = {
    'Spermatogonia': ['FGFR3', 'GAGE12J'],
    'Spermatocytes': ['DDX4', 'SYCP3',  'ZPBP'], # 'SPO11',
    'Spermatid': ['PRM3', 'LINC00467', 'TSACC'],
    'Sertoli': ['SOX9', 'MIR202HG',], # C4
    'Peritubular Myoid and Leydig': ['ACTA2', 'IGF1'],
    'Smooth Musc.': ['DCN'],
    'Mac': ['CD14', 'S100A4', 'FCER1G'],
    'Endothelial': ['VWF', 'GNG11', 'IFI27'],

    'markers' : ['UTF1', 'NODAL', 'GDF3', 'DPPA5', 'LEFTY1', 'LEFTY2', 'DPPA2',
    # HERVE-splice candidates:
    'SSUH2',
    'AC079466.1',
    'AC092666.1',
    'LINC01198',
    'SSUH2',
    'AP002495.1',
    ]

    }


# In[8]:


# ovary_marker3 (add ddx4)

marker_genes_dict = {
    #'Epiblast' [A1BG, POU5F1, 'NANOG', 'ERVB7_1-LTR_MM', 'ETnERV-int', 'RLTR13G'],
    #'Primitive streak' [Eomes, Nanog], #Nanong!!
    #'Anterior primitive streak' [Gsc, Mixl1],
    #'Notochord' [Noto, T],
    #'Def. Endoderm' [Cer1, Sox7],
    #'Nascent mesoderm' [Mesp1, Apela],
    #'Caudal mesoderm' [Cdx1, Hes7],
    #'Paraxial mesoderm' [Tcf15, Tbx1],
    #'Somitic mesoderm' [Tbx6, Dll1, 'RLTR1D2_MM'],
    #'Pharngyeal mesoderm' [Tcf21, Isl1],
    #'Cardiomyocytes' [Tnnt2, Myl4, 'RLTR13A2', 'ETnERV3-int', 'L1ME3D', 'L1MEc', 'RLTR16'],
    #'Allantois' [Tbx4, Hoxa11],
    #'Mesenchyme' [Krt18, Pmp22],
    #'Hemandothelial prog.' [Kdr, Etv2],
    #'Endothelium' [Pecam1, Anxa5],
    #'Blood prog.' [Runx1, Lmo2],
    #'Erythroid' [Gata1, Gypa, 'L1_Mur1', 'Lx2B2', 'RLTR10F'],
    #'Neuromesoderml prog.' [Cdx4, Epha5, 'ERVB4_1C-LTR_Mm', 'RMER17C2'],
    #'Neurectoderm' [Six3, Irx3],
    #'Neural crest' [Dlx2, Sox10],
    #'Brain' [En1, Pax2],
    #'Spinal cord' [Sox2, Pax2],
    #'Surface ectoderm' [Trp63, Grhl2],
    #'Visceral endoderm' [Dkk1, Amot],
    # 'Exe endoderm' [Ttr, Apoa2, 'RLTR1B-int', 'MER5C', 'Charlie10b', 'LTRIS2', 'MER46C', 'MuRRS4-int', 'LTR84b', 'RLTR20B3'],
    # 'Exe ectoderm' [Tfap2c, Elf5, 'RLTR45', 'RLTR9E', 'RLTR45-int', 'ERVB4_2-LTR_MM'],
    #'Parietal endoderm' [Sparc, Plat],
    #'Unknown' ['RLTR13C2', 'RMER17A-int',],
    'test': ["DDX6", "TFAP4", "MGA", "KIF13A", "AluY", "C1QTNF9B"]
    }


# In[9]:


sc.pl.stacked_violin(adata, marker_genes_dict, groupby='leiden_r0.10', rotation=90, dendrogram=True, show=False, save='markers.pdf')
sc.pl.dotplot(adata, marker_genes_dict, groupby='leiden_r0.10', dot_max=0.5, dendrogram=True, standard_scale='var', show=False, save='markers.pdf')
sc.pl.heatmap(adata, marker_genes_dict, groupby='leiden_r0.10', show=False, save='markers.pdf')

for k in marker_genes_dict:
    sc.pl.tsne(adata, color=marker_genes_dict[k], size=10, legend_loc='on data', vmax=3, show=False, save='markers-{0}.pdf'.format(k))
    sc.pl.umap(adata, color=marker_genes_dict[k], color_map='plasma', size=10, vmax=3, legend_loc='on data', show=False, save='markers-{0}.pdf'.format(k))


# In[ ]:




