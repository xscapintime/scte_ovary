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

sc.settings.figdir = 'markers'

sc.tl.rank_genes_groups(adata, 'leiden_r0.50', method='wilcoxon')
#sc.tl.filter_rank_genes_groups(adata, min_fold_change=2, min_in_group_fraction=0.25)
sc.pl.rank_genes_groups_dotplot(adata, key='rank_genes_groups', show=False, save='l0.5.pdf')

sc.tl.rank_genes_groups(adata, 'leiden_r0.20', method='wilcoxon')
#sc.tl.filter_rank_genes_groups(adata, min_fold_change=2, min_in_group_fraction=0.25)
sc.pl.rank_genes_groups_dotplot(adata, key='rank_genes_groups', show=False, save='l0.2.pdf')

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

# Selected From the powell paper:
markers = ['UTF1', 'NODAL', 'GDF3', 'DPPA5', 'LEFTY1', 'LEFTY2', 'DPPA2',
    # HERVE-splice candidates:
    'SSUH2',
    'AC079466.1',
    'AC092666.1',
    'LINC01198',
    'SSUH2',
    'AP002495.1',
    ]

print(adata)

sc.pl.umap(adata, color=markers, size=10, legend_loc='on data', show=False, save='-powell.pdf', vmin=0, vmax=2)
#sc.pl.stacked_violin(adata, var_names=gene, groupby='leiden_r0.2', rotation=90, show=False, save='{0}.pdf'.format(gene))
#print(gene)



