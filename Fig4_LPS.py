import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats

# figure params
sc.set_figure_params(figsize=(4, 4), fontsize=15)

# reading in preprocessing data object, all four conditions
bmdm_all = sc.read('/Users/katebridges/Downloads/bmdm_object_raw.h5ad')

# removing m1+m2 condition for fig 3
bmdm_fig3 = bmdm_all[np.where(bmdm_all.obs['sample'] != 'M1+M2')[0], :]

# 2D embedding of full dataset
del bmdm_fig3.obsm, bmdm_fig3.obsp, bmdm_fig3.uns

# working toward visualization
sc.tl.pca(bmdm_fig3, svd_solver='auto')
sc.pl.pca(bmdm_fig3, color='sample')
sc.pp.neighbors(bmdm_fig3)  # using with default parameters

# UMAP VIZ of three expr conditions
sc.tl.umap(bmdm_fig3)  # default initial position
sc.pl.umap(bmdm_fig3, color='sample')
sc.tl.leiden(bmdm_fig3, resolution=0.2)

# inspecting expression of GENES OF INTEREST
kmj_genes = ['Il12b', 'Il6', 'Tnf', 'Ccl5', 'H2-Eb1', 'H2-Aa', 'H2-Ab1', 'Ccr7', 'Flt3', 'Ccl22', 'Ccl17', 'Socs2']
for j in kmj_genes:
    sc.pl.umap(bmdm_fig3, color=j, cmap='Reds')

# PLOTTING bar graphs of genes oi w bootstrapped error bars
clust2cond = {'0': 'BMDM Control',
              '1': 'BMDM IL-4',
              '2': 'BMDM LPS+IFNg',
              '3': 'BMDM LPS+IFNg',
              '4': 'BMDM Control',
              '5': 'Cluster 5'}

clust_lab = []
for j in bmdm_fig3.obs['leiden'].values:
    clust_lab.append(clust2cond[j])

bmdm_fig3.obs['clust2cond'] = clust_lab

bar_ind = np.unique(bmdm_fig3.obs['clust2cond'])
for gene_name in kmj_genes:
    dat = np.array(bmdm_fig3[:, gene_name].X.todense()).flatten()
    dat_stat = np.zeros((len(np.unique(bmdm_fig3.obs['clust2cond'])), 3))
    b = 0
    for g in np.unique(bmdm_fig3.obs['clust2cond']):
        i = np.where(bmdm_fig3.obs['clust2cond'] == g)[0]
        ci_info = bs.bootstrap(dat[i], stat_func=bs_stats.mean)
        dat_stat[b, 0] = ci_info.value
        dat_stat[b, 1] = dat_stat[b, 0] - ci_info.lower_bound
        dat_stat[b, 2] = ci_info.upper_bound - dat_stat[b, 0]
        b = b + 1
    fig, ax = plt.subplots(figsize=(4, 5.5))
    barlist = ax.bar(bar_ind, dat_stat[:, 0], yerr=[dat_stat[:, 1], dat_stat[:, 2]], align='center', ecolor='black', capsize=10)
    barlist[0].set_color(sns.husl_palette(6)[4])
    barlist[1].set_color(sns.husl_palette(6)[2])
    barlist[2].set_color('C1')
    barlist[3].set_color(sns.husl_palette(6)[5])
    plt.title(gene_name)
    ax.set_ylabel('ln[mRNA counts + 1]')
    ax.set_xticks(np.arange(len(bar_ind)))
    ax.set_xticklabels(bar_ind, rotation=60)
    # ax.set_xlim([-0.5, len(bar_ind) + 0.5])
    plt.tight_layout()

# DIFF EXPR analysis for cluster 5 vs. rest after LPS stim
# distinguishing cluster 5 from rest in metadata
sample_cluster = []
for k in range(bmdm_fig3.shape[0]):
    if bmdm_fig3.obs['leiden'][k] == '5':
        sample_cluster.append(bmdm_fig3.obs['sample'][k] + ' ' + bmdm_fig3.obs['leiden'][k])
    else:
        sample_cluster.append(bmdm_fig3.obs['sample'][k])

bmdm_fig3.obs['sample_cluster'] = sample_cluster

# limiting evaluation to M1 condition only
bmdm_m1 = bmdm_fig3[np.where(bmdm_fig3.obs['sample'] == 'M1')[0], :]

# calculating diff expr genes
sc.tl.rank_genes_groups(bmdm_fig3, groupby='sample_cluster', groups=['M1 5', 'M1'],
                        key_added='LPS')
dc_lps = sc.get.rank_genes_groups_df(bmdm_fig3, group='M1 5', key='LPS', pval_cutoff=0.05, log2fc_min=1.5)
mac_lps = sc.get.rank_genes_groups_df(bmdm_fig3, group='M1', key='LPS', pval_cutoff=0.05, log2fc_min=1.5)

# write results to xlsx -> for import to gProfiler for GSEA
excelpath = '/Users/katebridges/PyCharmProjects/test/bmdm_lps_cluster5_20211025.xlsx'
writer = pd.ExcelWriter(excelpath, engine='xlsxwriter')

dc_lps.to_excel(writer, sheet_name='Cluster5+LPS')
mac_lps.to_excel(writer, sheet_name='Rest+LPS')

writer.save()
