# %%
import logging, matplotlib, os, sys
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors

# %%
from glbase3 import *

# %%
plt.rcParams['figure.figsize']=(8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['font.size']=6

# %%
sc.settings.figdir = 'specific-tes'

# %%
#%%time

adata = sc.read('./learned.h5ad')

# %%
marker_genes_dict = {
    'test':["FIGLA"]
    #'Epiblast': ["Pou5f1", 'Nanog', 'RodERV21-int', 'ERVB7_1-LTR_MM', 'ETnERV-int', 'RLTR13G'],
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
    #'Exe endoderm': ["Ttr", "Apoa2", 'RLTR1B-int', 'MER5C', 'Charlie10b', 'LTRIS2', 'MER46C', 'MuRRS4-int', 'LTR84b', 'RLTR20B3'],
    #'Exe ectoderm': ["Tfap2c", "Elf5", 'RLTR45', 'RLTR9E', 'RLTR45-int', 'ERVB4_2-LTR_MM'],
    #'Parietal endoderm': ["Sparc", "Plat"],
    #'Unknown': ['RLTR13C2', 'RMER17A-int',],
    }

# %%
# need find the specific marker genes
# the code works

# marker_genes_dict = {
#      'test': ['AAAS'],
# }

# %%
#%%time

sc.pl.stacked_violin(adata, marker_genes_dict, groupby='leiden_r0.50', rotation=90, dendrogram=True, show=False, save='markers.pdf')
sc.pl.dotplot(adata, marker_genes_dict, groupby='leiden_r0.50', dot_max=0.7, dendrogram=True, standard_scale='var', vmax=2, show=False, save='markers.pdf')
sc.pl.heatmap(adata, marker_genes_dict, groupby='leiden_r0.50', show=False, save='markers.pdf', vmax=3)

# %%
sc.pl.stacked_violin(adata, marker_genes_dict, groupby='leiden_r0.20', rotation=90, dendrogram=True, show=False, save='markers-0.2.pdf')
sc.pl.dotplot(adata, marker_genes_dict, groupby='leiden_r0.20', dot_max=0.7, dendrogram=True, standard_scale='var', vmax=1, show=False, save='markers-0.2.pdf')
sc.pl.heatmap(adata, marker_genes_dict, groupby='leiden_r0.20', show=False, save='markers-0.2.pdf', vmax=3)

# %%
for k in marker_genes_dict:
    sc.pl.tsne(adata, color=marker_genes_dict[k], size=10, legend_loc='on data', vmax=2, show=False, save='markers-{0}.pdf'.format(k))
    sc.pl.umap(adata, color=marker_genes_dict[k], size=10, legend_loc='on data', vmax=2, show=False, save='markers-{0}.pdf'.format(k))

# %%
print(adata.obs)

# %%
print(adata.raw.var_names)

# %%



