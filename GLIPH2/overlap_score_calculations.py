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
#%% load gliph results *GLIPH*

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
tissue_high_grouping = {'large_intestine': 'GI_track', 'small_intestine': 'GI_track',  'stomach': 'GI_track',
                        'esophagus': 'GI_track', 'skin': 'skin', 'heart': 'control', 'kidney': 'control',
                        'liver': 'liver', 'blood': 'blood', 'LN': 'LN', 'spleen' : 'spleen', 'pretx_blood': 'pretx_blood'}
#%%
gliph_all = pd.read_csv('gliph/test_cluster_allpts_withcomparators.csv')
gliph_all[['Patient','Tissue']] = gliph_all.Sample.str.split(':',expand=True)
gliph_all.rename(columns={'index': 'gliph_index'}, inplace=True)

gliph_all['Tissue_group'] = gliph_all['Tissue'].map(tissue_grouping)
gliph_all['Tissue_high_group'] = gliph_all['Tissue_group'].map(tissue_high_grouping)



#%%

conds = [(gliph_all[' number_unique_cdr3']>2,'number_unique_cdr3>2'),
         (gliph_all.Fisher_score<0.05, 'Fisher<0.05'),
         (gliph_all.vb_score<0.05, 'vb_score<0.05'),
         (gliph_all.length_score<0.05, 'length_score<0.05'),
         ]

#%% filtering gliph table
cond = conds[0:1] #choose which condition to use for filtering
gliph_filter = gliph_all[np.logical_and.reduce([x[0] for x in cond])]
T = ' '.join([x[1] for x in cond]) + ' n=' + str(len(gliph_filter.gliph_index.unique()))

#%% filtering gliph table

#gliph_filter = gliph_all[(gliph_all.Fisher_score<0.05) & (gliph_all.vb_score<0.05) & (gliph_all.length_score<0.05) & (gliph_all[' number_unique_cdr3']>2)]
#T = 'Fisher<0.05 vb_score<0.05 length_score<0.05 number_unique_cdr3>2 n=' + str(len(gliph_filter.gliph_index.unique()))


#%% list of gliph groups

gliph_clusters = pd.concat([gliph_filter.iloc[:,np.r_[:12]].groupby('gliph_index').agg('first'),
                            gliph_filter.iloc[:,np.r_[0,29:33]].groupby('gliph_index').agg('nunique'),
                            gliph_filter.iloc[:,0:1].groupby('gliph_index').size().rename('cluster_size')], axis=1)

#%%
iix_n = pd.MultiIndex.from_product([tissue_groups_order, np.unique(gliph_all.Patient)])
A = gliph_filter.pivot_table(values='pattern', index='gliph_index', columns=['Tissue_group','Patient'],
                          aggfunc=pd.Series.nunique).reindex(columns=iix_n)

# boolean 3d matrix gliph/tissue/patient for every combination that exist
combs = A.fillna(0).to_numpy().reshape(gliph_filter.gliph_index.nunique(),
                                       gliph_filter.Tissue_group.nunique(),
                                       gliph_filter.Patient.nunique())
# summing over tissue for total sum for each group in each patient
B = A.groupby(level=1, axis=1).sum()

# sum over patients - how many patients have each gliph/tissue comb
C = np.sum(combs, axis=2)
C_df = pd.DataFrame(data=C, columns=tissue_groups_order, index=A.index, dtype=int)

#%%

gliph_clusters['GVHD'] = B[['Pt1', 'Pt2', 'Pt3', 'Pt5', 'Pt6', 'Pt7', 'Pt9']].sum(axis=1)>0
gliph_clusters['comp'] = B[['Pt10','Pt11','Pt12']].sum(axis=1)>0

gliph_clusters['GVHD_only'] = gliph_clusters['GVHD'] & ~gliph_clusters['comp']
gliph_clusters['comp_only'] = gliph_clusters['comp'] & ~gliph_clusters['GVHD']
gliph_clusters['GVHD_comp'] = gliph_clusters['comp'] & gliph_clusters['GVHD']

gliph_clusters['patients'] = ''
gliph_clusters['patients'][gliph_clusters['GVHD_only']] = 'GVHD only'
gliph_clusters['patients'][gliph_clusters['comp_only']] = 'no GVHD'
gliph_clusters['patients'][gliph_clusters['GVHD_comp']] = 'mixed'

gliph_clusters['length'] = gliph_clusters['pattern'].str.len()
gliph_clusters['type_short'] = gliph_clusters['type'].str.split('-').str[0]
#%% interesting gliph groups
interest_size = [6, 3, 3, 4, 4, 7, 6, 3, 3, 2, 2, 2]
interesting_groups = (C>=interest_size).any(axis=1)


pd.concat([gliph_clusters[interesting_groups], C_df[interesting_groups]],axis=1).to_csv('interesting_gliph_groups_filtered.csv')



#%% calculating sample sizes and aggragated frequenices (sum of freq for each gliph group and sample)

sample_sizes2 = gliph_all[['TcRb','V','J','Freq','Patient',
                          'Tissue_group', 'Sample']].drop_duplicates().groupby(['Tissue_group','Patient'])['Freq'].sum()


agg_freq = gliph_filter.pivot_table(values='Freq', index='gliph_index', columns=['Tissue_group','Patient'],
                                 aggfunc=pd.Series.sum).fillna(0) / sample_sizes2
#%% calculating Jaccard overlap
mapper_dict = {x: int(x[2:]) for x in agg_freq.columns.levels[1]}
jacc = pd.DataFrame(1 - squareform(pdist(agg_freq.T.astype(bool), metric='jaccard')),
                        index=agg_freq.columns,
                        columns=agg_freq.columns).rename(index=mapper_dict).rename(columns=mapper_dict).sort_index(level=0).sort_index(axis=1, level=0)



#%%
tissue_to_plot = ['A_blood', 'spleen', 'esophagus', 'stomach', 'small_intestine', 'large_intestine', 'liver', 'skin']
#groups_to_plot = very_interesting_groups#agg_freq.loc[:,tissue_to_plot].max(axis=1) > 0.1
patients_to_plot = [1, 2, 3, 5, 6, 7, 9, 10, 11, 12]

jacc_to_plot = jacc.loc[(tissue_to_plot,patients_to_plot),(tissue_to_plot,patients_to_plot)].rename(
    {'A_blood': 'blood', 'large_intestine': 'large intestine', 'small_intestine': 'small intestine'}).rename(
    {'A_blood': 'blood', 'large_intestine': 'large intestine', 'small_intestine': 'small intestine'}, axis=1)

long_jacc_to_plot = jacc_to_plot.reset_index().rename(
    columns={'Tissue_group':'Tissue_group_1', 'Patient' : 'Patient_1'}).melt(
    id_vars=['Tissue_group_1', 'Patient_1']).rename(
    columns={'Tissue_group':'Tissue_group_2', 'Patient' : 'Patient_2'})
#%% mask same patient samples
M = pd.DataFrame(index=jacc_to_plot.index, columns=jacc_to_plot.index)

for i in M.index:
    for j in M.index:
        M.loc[i,j] = True if i[1]==j[1] else False

#%% heatmap of gliph jaccard by tissue
plt.rcParams.update({'font.size': 10})
plt.figure(figsize=(16,16))
n = len(jacc_to_plot)
ax = sns.heatmap(jacc_to_plot, mask=np.triu(np.ones([n,n])).astype(bool) | M, cmap='Blues',#center=0.04 ,cmap='PRGn',
                 square=True, cbar_kws={"shrink": .6})
ax.set_facecolor('lightgrey') #lightgrey
for t in np.cumsum(jacc_to_plot.groupby(level=0).nunique().iloc[:,0].reindex(tissue_to_plot).to_list()):
    ax.axhline(y=t, c='white')
    ax.axvline(x=t, c='white')
plt.tight_layout()
plt.savefig('gliph_jacard_by_tis.pdf')
#plt.show()

#%% heatmap of gliph jaccard by patient
plt.figure(figsize=(16,16))
n = len(jacc_to_plot)
ax = sns.heatmap(jacc_to_plot.sort_index(level=1).sort_index(axis=1, level=1),
            mask=np.triu(np.ones([n,n])).astype(bool), vmax=0.3, center=0.1, cmap='Blues')
ax.set_facecolor('lightgrey')
for t in np.cumsum(jacc_to_plot.groupby(level=1).nunique().iloc[:,0].to_list()):
    ax.axhline(y=t, c='lightgrey')
    ax.axvline(x=t, c='lightgrey')
plt.tight_layout()
plt.savefig('gliph_jacard_by_pt.pdf')
#%%
sample_sizes_ptnum = sample_sizes2.rename(index=mapper_dict).sort_index()
sample_sizes_prod = pd.DataFrame(np.outer(sample_sizes_ptnum, sample_sizes_ptnum),
                                 index=sample_sizes_ptnum.index, columns=sample_sizes_ptnum.index)
sample_sizes_geomean = np.sqrt(sample_sizes_prod)

#%%
J_long = jacc.where(np.triu(
    np.ones(jacc.shape), 1).astype(bool)).melt(
    ignore_index=False, value_name='jaccard').rename(
    columns={'Tissue_group' : 'T2', 'Patient' : 'P2'}).reset_index().rename(
    columns={'Tissue_group' : 'T1', 'Patient' : 'P1'}).dropna()

S_long = sample_sizes_geomean.where(np.triu(
    np.ones(sample_sizes_geomean.shape), 1).astype(bool)).melt(
    ignore_index=False, value_name='sample_size').rename(
    columns={'Tissue_group' : 'T2', 'Patient' : 'P2'}).reset_index().rename(
    columns={'Tissue_group' : 'T1', 'Patient' : 'P1'}).dropna()

JS_long = pd.concat([J_long, S_long], axis=1).loc[:,~pd.concat([J_long, S_long], axis=1).columns.duplicated()]

JS_long['num_comp'] = (JS_long.P1>9).astype(int) + (JS_long.P2>9).astype(int)
JS_long['num_pt5'] = (JS_long.P1==5).astype(int) + (JS_long.P2==5).astype(int)
JS_long['same_patient'] = (JS_long.P1==JS_long.P2)
JS_long['num_int'] = (JS_long.T1.str.endswith('tine')).astype(int) + (JS_long.T2.str.endswith('tine')).astype(int)
GI = ['stomach', 'large_intestine', 'small_intestine', 'esophagus']
JS_long['num_GI'] = (JS_long.T1.isin(GI)).astype(int) + (JS_long.T2.isin(GI)).astype(int)

JS_long['Ts'] = JS_long[['T1', 'T2']].values.tolist()
JS_long['tissues'] = JS_long.T1 + '/' +  JS_long.T2

#%%
JS_long['jaccard_log'] = np.log10(JS_long['jaccard'])
JS_long['sample_size_log'] = np.log10(JS_long['sample_size'])
JS_long.jaccard_log[(JS_long.jaccard==0)] = None

#%% jaccard vs sample size comparators
#plt.figure(figsize=(12,8))
f = sns.lmplot(data=JS_long[(JS_long.T1 != 'pretx_blood') & (JS_long.T2 != 'pretx_blood') & (JS_long.num_pt5 == 0) & (JS_long.same_patient == 0)],
                x='sample_size_log', y='jaccard_log', hue='num_comp', truncate=False, height=8, aspect=1.2)
plt.xlabel('log of geometric mean of sample sizes')
plt.ylabel('log of Jaccard similarity of Gliph groups')
#plt.legend(title='Number of comparators in pair')
f._legend.set_title('Number of \ncomparators in pair')
plt.savefig('comparators_jaccard_reg.pdf')

#%% jaccard vs sample size pt5
#plt.figure(figsize=(12,8))
f = sns.lmplot(data=JS_long[(JS_long.T1 != 'pretx_blood') & (JS_long.T2 != 'pretx_blood') & (JS_long.same_patient == 0)],
                x='sample_size', y='jaccard', hue='num_pt5', palette=['red','blue'], height=8, aspect=1.2, fit_reg=False)#, logx=True)
f.set(xscale="log")
f.set(yscale="log")
#plt.ylim([1e-3,1e0])
#plt.xlim([1e4,2e10])
plt.savefig('pt5_jaccard_reg.pdf')
#%%
plt.figure(figsize=(12,8))
f = sns.relplot(data=JS_long, x='sample_size', y='jaccard', hue='num_pt5', col='same_patient')
f.set(xscale="log")
f.set(yscale="log")
plt.suptitle(T)
#plt.ylim([1e-3,1e0])
#plt.xlim([1e4,2e10])

#plt.show()

#%% universal slope of jaccard vs sample size
d = JS_long[(JS_long.T1 != 'pretx_blood') & (JS_long.T2 != 'pretx_blood')
            & (JS_long.num_pt5 == 0) & (JS_long.same_patient == 0)]
x = d['sample_size_log']
y = d['jaccard_log']

x = x[~y.isna()]
y = y[~y.isna()]

fit, covar = curve_fit(lambda x,a,b : a*x + b, x, y)
slope = fit[0]

# %% standard sharing
x_stan = 3

#%% calculating standard interception

intercept_all = JS_long[~JS_long['Ts'].str.contains('pretx_blood', regex=False) &
                        (JS_long.num_pt5 == 0) &
                        (JS_long.same_patient == 0)].groupby(['T1', 'T2']).apply(lambda d :
                                                              curve_fit(lambda x,b :
                                                                        slope * (x-x_stan) + b,
                                                                        d['sample_size_log'][~d['jaccard_log'].isna()],
                                                                        d['jaccard_log'][~d['jaccard_log'].isna()])[0][0]).unstack(level=1)

intercept_sym = intercept_all.fillna(0) + intercept_all.fillna(0).T - \
                pd.concat([intercept_all.fillna(0), intercept_all.fillna(0).T]).max(level=0)#%%

intercept_sym_zero_diag = (abs(intercept_sym) - np.min(np.array(abs(intercept_sym))))/(np.max(np.array(abs(intercept_sym))) - np.min(np.array(abs(intercept_sym))))
intercept_sym_zero_diag = intercept_sym_zero_diag - np.diag(np.diag(intercept_sym_zero_diag))



#%%
Z = linkage(squareform(intercept_sym_zero_diag))
#plt.figure(figsize=(8, 8))
f = sns.clustermap(intercept_sym,
               annot=True, fmt=".3", row_linkage=Z, col_linkage=Z, cmap='PuOr_r') #'seismic_r'
plt.title(T)
#f.savefig('interecepts_clustering_full_'+ T +'.pdf')

mask = np.triu(np.ones([8, 8]),k=1).astype(bool)
mask = mask[np.argsort(f.dendrogram_row.reordered_ind),:]
mask = mask[:,np.argsort(f.dendrogram_col.reordered_ind)]

f = sns.clustermap(intercept_sym,
               annot=True, linewidth=.5, fmt=".3", row_linkage=Z, col_linkage=Z, cmap='PuOr_r',
                   mask=mask)
plt.tight_layout()
plt.title(T)
f.savefig('interecepts_clustering_tri_'+ T +'.pdf')

#%% plotting pairs with one SI with constrained slope fit
plt.figure(figsize=(12,8))
f = sns.scatterplot(data=JS_long[~JS_long['Ts'].str.contains('pretx_blood', regex=False) &
                            JS_long['Ts'].str.contains('small_intestine', regex=False) &
                            (JS_long.num_pt5 == 0) & (JS_long.same_patient == 0)],
                x='sample_size_log', y='jaccard_log', hue='tissues', legend=False)
sns.lineplot(data=pd.DataFrame(data=np.add.outer(slope*(np.array(f.get_xlim())-x_stan),
                                            np.array(intercept_sym['small_intestine'])),
                          index=f.get_xlim(), columns=intercept_sym.index).rename(
    columns = {'A_blood': 'blood', 'large_intestine': 'large intestine', 'small_intestine': 'small intestine'}).stack().reset_index().rename(
    columns = {'level_0': 'x', 0: 'yy'}), x='x', y='yy', hue='T1')

plt.legend(title = 'tissues')
#plt.legend(loc='lower right', title='tissues')
plt.savefig('reg_SI_fixed_slope_vs_rest.pdf')
