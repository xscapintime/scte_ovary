#!/usr/bin/env python
# coding: utf-8

# In[2]:


# %load 8.detect_markers.py
import logging, matplotlib, os, sys
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
plt.rcParams['figure.figsize']=(8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['font.size']=10
adata = sc.read('./learned.h5ad')


# In[3]:


sc.settings.figdir = 'markers'

sc.tl.rank_genes_groups(adata, 'leiden_r0.50', method='wilcoxon')
#sc.tl.filter_rank_genes_groups(adata, min_fold_change=2, min_in_group_fraction=0.25)
sc.pl.rank_genes_groups_dotplot(adata, key='rank_genes_groups', show=False, save='l0.5.pdf')

sc.tl.rank_genes_groups(adata, 'leiden_r0.20', method='wilcoxon')
#sc.tl.filter_rank_genes_groups(adata, min_fold_change=2, min_in_group_fraction=0.25)
sc.pl.rank_genes_groups_dotplot(adata, key='rank_genes_groups', show=False, save='l0.2.pdf')


# In[4]:


p = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
print(p.head(5))
p.to_csv('marker_genes_per_cluster.tsv', sep='\t')

parts = {1: slice(0, 10), 2: slice(10, 20), 3: slice(20, 30), 4: slice(30, 40)}

for grp in p.columns.values:
    print(grp)
    col = p[grp]

    for k in parts:
        sc.pl.umap(adata, color=col[parts[k]], size=10, legend_loc='on data', show=False, save='-g{0}-p{1}.pdf'.format(grp, k), vmin=0, vmax=2)
        sc.pl.dotplot(adata, col[parts[k]], groupby='leiden_r0.20', show=False, save='-g{0}-p{1}.pdf'.format(grp, k))


# In[12]:


# Selected From the powell paper:
# delete 
markers = ['UTF1', 'NODAL', 'GDF3', 'DPPA5', 'LEFTY1', 'LEFTY2', 'DPPA2',
    # HERVE-splice candidates:
    'SSUH2',
    'AC079466.1',
    'AC092666.1',
    'LINC01198',
    'AP002495.1',
    ]

print(adata)


# In[20]:


# ovary_marker3 (add ddx4)

markers = {
    #'Epiblast': ["A1BG", "POU5F1", 'NANOG', 'ERVB7_1-LTR_MM', 'ETnERV-int', 'RLTR13G'],
    #'Primitive streak': ["Eomes", "Nanog"], #Nanong?!?!
    #'Anterior primitive streak': ["Gsc", "Mixl1"],
    #'Notochord': ["Noto", "T"],
    #'Def. Endoderm': ["Cer1", "Sox7"],
    #'Nascent mesoderm': ["Mesp1", "Apela"],
    #'Caudal mesoderm': ["Cdx1", "Hes7"],
    #'Paraxial mesoderm': ["Tcf15", "Tbx1"],
    #'Somitic mesoderm': ["Tbx6", "Dll1", 'RLTR1D2_MM'],
    #'Pharngyeal mesoderm': ["Tcf21", "Isl1"],
    #'Cardiomyocytes': ["Tnnt2", "Myl4", 'RLTR13A2', 'ETnERV3-int', 'L1ME3D', 'L1MEc', 'RLTR16'],
    #'Allantois': ["Tbx4", "Hoxa11"],
    #'Mesenchyme': ["Krt18", "Pmp22"],
    #'Hemandothelial prog.': ["Kdr", "Etv2"],
    #'Endothelium': ["Pecam1", "Anxa5"],
    #'Blood prog.': ["Runx1", "Lmo2"],
    #'Erythroid': ["Gata1", "Gypa", 'L1_Mur1', 'Lx2B2', 'RLTR10F'],
    #'Neuromesoderml prog.': ["Cdx4", "Epha5", 'ERVB4_1C-LTR_Mm', 'RMER17C2'],
    #'Neurectoderm': ["Six3", "Irx3"],
    #'Neural crest': ["Dlx2", "Sox10"],
    #'Brain': ["En1", "Pax2"],
    #'Spinal cord': ["Sox2", "Pax2"],
    #'Surface ectoderm': ["Trp63", "Grhl2"],
    #'Visceral endoderm': ["Dkk1", "Amot"],
    # 'Exe endoderm': ["Ttr", "Apoa2", 'RLTR1B-int', 'MER5C', 'Charlie10b', 'LTRIS2', 'MER46C', 'MuRRS4-int', 'LTR84b', 'RLTR20B3'],
    # 'Exe ectoderm': ["Tfap2c", "Elf5", 'RLTR45', 'RLTR9E', 'RLTR45-int', 'ERVB4_2-LTR_MM'],
    #'Parietal endoderm': ["Sparc", "Plat"],
    #'Unknown': ['RLTR13C2', 'RMER17A-int',],
    "DDX6", "TFAP4", "MGA", "KIF13A", "AluY", "C1QTNF9B"
    }


# In[21]:


sc.pl.umap(adata, color=markers, size=10, legend_loc='on data', show=False, save='-powell.pdf', vmin=0, vmax=2)
#sc.pl.stacked_violin(adata, var_names=gene, groupby='leiden_r0.2', rotation=90, show=False, save='{0}.pdf'.format(gene))
#print(gene)

