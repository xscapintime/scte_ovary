import os, sys, glob
import pandas as pd



filtred_bc = pd.read_csv('adata_srr.csv', header=0, sep=',')
# filtred_bc['bc_srr'] = filtred_bc['barcodes'] + '_' + filtred_bc['srr']
sel_sample = filtred_bc['srr'].unique().tolist()

# srr_sampl = {sam for sam in filtred_bc['sample_id'].unique().tolist()}
# srr_sampl = pd.read_csv('../fan_project/metadata/metadata.txt', header=None, sep='\t')
# sampl_dict = dict(srr_sampl.values)
# filtred_bc['srr'] = filtred_bc['sample_id'].map(sampl_dict)


path = '../fan_project/fan2019_data'
files = glob.glob(os.path.join(path, '*.gz'))
files =  [x for x in files if os.path.split(x)[1].split('_')[0] in sel_sample ]



df_ls = []

for file in files[8:10]:
    df = pd.read_csv(file, header=0, sep=',', compression='gzip')
    df['srr'] = os.path.split(file)[1].split('_')[0]
    # df['bc_srr'] = df['barcodes'] + '_' + df['srr']

    ## after filtering, some df will have 0 rows left
    # merged = pd.merge(df, filtred_bc, on=['barcodes', 'srr'], how='inner')        
    # trans = merged.drop(labels=['cluster', 'srr'], axis=1).set_index('barcodes').transpose()
    
    # names = os.path.split(file)[1].split('_')[0]
    df_ls.append(df)



all_sample = pd.concat(df_ls, axis=0)
filtered = pd.merge(all_sample, filtred_bc, on=['barcodes', 'srr'], how='inner')
filtered['idx'] = filtered['barcodes'] + '_' + filtered['srr'] + '_' + filtered['cluster']
fin = filtered.drop(labels=['cluster', 'srr', 'barcodes'], axis=1).set_index('idx').transpose()


fin.to_csv('allsam_T.txt.gz', sep='\t', header=True, index=True, compression='gzip')