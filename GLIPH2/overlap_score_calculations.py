import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib.pylab as plt
import matplotlib
import sys
from scipy.stats import entropy
from scipy.optimize import curve_fit
#%%
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage

#%% load gliph results



tissue_grouping = {'A_blood': 'A_blood', 'A_marrow': 'A_blood', 'ascending_colon': 'large_intestine',
                    'descending_colon': 'large_intestine', 'ileum': 'small_intestine', 'spleen' : 'spleen',
                    'liver': 'liver', 'jejunum': 'small_intestine', 'rectum': 'large_intestine',
                    'skin': 'skin', 'stomach': 'stomach', 'duodenum': 'small_intestine',
                    'esophagus': 'esophagus', 'transverse_colon': 'large_intestine', 'pretx_blood': 'pretx_blood',
                    'left_colon': 'large_intestine', 'right_colon': 'large_intestine',
                    'colon': 'large_intestine', 'L_colon': 'large_intestine', 'kidney': 'kidney',
                    'mLN': 'mLN', 'mid_colon': 'large_intestine', 'small_intestine': 'small_intestine',
                    'heart': 'heart', '2_skin': 'skin'}

tissue_groups_order = ['pretx_blood','A_blood', 'spleen', 'esophagus', 'stomach', 'small_intestine', 'large_intestine',
                       'liver', 'skin', 'heart', 'kidney', 'mLN']

#%%
gliph_all = pd.read_csv(sys.argv[1])
gliph_all[['Patient','Tissue']] = gliph_all.Sample.str.split(':',expand=True)
gliph_all.rename(columns={'index': 'gliph_index'}, inplace=True)

gliph_all['Tissue_group'] = gliph_all['Tissue'].map(tissue_grouping)

#%%
# calcuates total number of TCRs in each group from gliph table
sample_sizes = gliph_all[['TcRb','V','J','Freq','Patient',
                          'Tissue_group', 'Sample']].drop_duplicates().groupby(['Tissue_group','Patient'])['Freq'].sum()
# uses pivot to create a table of TCR frequencies by patient, tissue type and gliph group
agg_freq = gliph_all.pivot_table(values='Freq', index='gliph_index', columns=['Tissue_group','Patient'],
                                 aggfunc=pd.Series.sum).fillna(0) / sample_sizes

mapper_dict = {x: int(x[2:]) for x in agg_freq.columns.levels[1]}
# matrix for the jaccard index based on gliph groups between each two patient/tissue combinations
jacc = pd.DataFrame(1 - squareform(pdist(agg_freq.T.astype(bool), metric='jaccard')),
                        index=agg_freq.columns,
                        columns=agg_freq.columns).rename(index=mapper_dict).rename(columns=mapper_dict).sort_index(level=0).sort_index(axis=1, level=0)

#%%
sample_sizes_ptnum = sample_sizes.rename(index=mapper_dict).sort_index()
sample_sizes_prod = pd.DataFrame(np.outer(sample_sizes_ptnum, sample_sizes_ptnum),
                                 index=sample_sizes_ptnum.index, columns=sample_sizes_ptnum.index)
sample_sizes_geomean = np.sqrt(sample_sizes_prod)  # geometric mean of sample sizes between pooled samples
#%%
#long format of Jacard index
J_long = jacc.where(np.triu(
    np.ones(jacc.shape), 1).astype(np.bool)).melt(
    ignore_index=False, value_name='jaccard').rename(
    columns={'Tissue_group' : 'T2', 'Patient' : 'P2'}).reset_index().rename(
    columns={'Tissue_group' : 'T1', 'Patient' : 'P1'}).dropna()
#long format of pooled sample sizes
S_long = sample_sizes_geomean.where(np.triu(
    np.ones(sample_sizes_geomean.shape), 1).astype(np.bool)).melt(
    ignore_index=False, value_name='sample_size').rename(
    columns={'Tissue_group' : 'T2', 'Patient' : 'P2'}).reset_index().rename(
    columns={'Tissue_group' : 'T1', 'Patient' : 'P1'}).dropna()

JS_long = pd.concat([J_long, S_long], axis=1).loc[:,~pd.concat([J_long, S_long], axis=1).columns.duplicated()]

# combined table of jaccard index and sample size for each pair of pooled samples
JS_long['num_pt5'] = (JS_long.P1==5).astype(int) + (JS_long.P2==5).astype(int)
JS_long['same_patient'] = (JS_long.P1==JS_long.P2)
JS_long['Ts'] = JS_long[['T1', 'T2']].values.tolist() # combained tissue for both samples
JS_long['jaccard_log'] = np.log10(JS_long['jaccard'])
JS_long['sample_size_log'] = np.log10(JS_long['sample_size'])
JS_long.jaccard_log[(JS_long.jaccard == 0)] = None

#%% calculating universal slope of sharing vs mean sample size
d = JS_long[(JS_long.T1 != 'pretx_blood') & (JS_long.T2 != 'pretx_blood')
            & (JS_long.num_pt5 == 0) & (JS_long.same_patient == 0)]
x = d['sample_size_log']
y = d['jaccard_log']

x = x[~y.isna()]
y = y[~y.isna()]

fit, covar = curve_fit(lambda x,a,b : a*x + b, x, y)
slope = fit[0]
#%% calculates standart sharing where the fit line meet 1000 mean sample size
x_stan = 3
intercept_all = JS_long[~JS_long['Ts'].str.contains('pretx_blood', regex=False) &
                        (JS_long.num_pt5 == 0) &
                        (JS_long.same_patient == 0)].groupby(['T1', 'T2']).apply(lambda d :
                                                              curve_fit(lambda x,b :
                                                                        slope * (x-x_stan) + b,
                                                                        d['sample_size_log'][~d['jaccard_log'].isna()],
                                                                        d['jaccard_log'][~d['jaccard_log'].isna()])[0][0]).unstack(level=1)

intercept_sym = intercept_all.fillna(0) + intercept_all.fillna(0).T - \
                pd.concat([intercept_all.fillna(0), intercept_all.fillna(0).T]).max(level=0)

#%% scaled sharing for plotting purposes

intercept_sym_zero_diag = (abs(intercept_sym) - 2.18)/0.31
intercept_sym_zero_diag = intercept_sym_zero_diag - np.diag(np.diag(intercept_sym_zero_diag))

#%%
Z = linkage(squareform(intercept_sym_zero_diag))
plt.figure(figsize=(8, 8))
f = sns.clustermap(intercept_sym,
               annot=True, fmt=".3", row_linkage=Z, col_linkage=Z, cmap='seismic_r')
mask = np.triu(np.ones([8, 8]),k=1).astype(bool)
mask = mask[np.argsort(f.dendrogram_row.reordered_ind),:]
mask = mask[:,np.argsort(f.dendrogram_col.reordered_ind)]

f = sns.clustermap(intercept_sym,
               annot=False, fmt=".3", row_linkage=Z, col_linkage=Z, cmap='seismic_r',
                   mask=mask)
plt.tight_layout()
plt.savefig("gliph_overlap.png")
