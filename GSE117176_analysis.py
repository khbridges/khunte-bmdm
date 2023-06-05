import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats

# figure params
sc.set_figure_params(figsize=(4, 4), fontsize=15)

bmdm = sc.read('/Users/katebridges/Downloads/bmdm_GSE117176.h5ad')

# towards visualization
sc.tl.pca(bmdm, svd_solver='auto')
sc.pl.pca(bmdm, color='sample')
sc.pp.neighbors(bmdm)  # using with default parameters

# UMAP visualization
sc.tl.umap(bmdm)  # default initial position
sc.tl.leiden(bmdm, resolution=0.2)

# need to compute DEGs for each leiden cluster & write to file for reference
# get DEGs for cluster of all data
sc.tl.rank_genes_groups(bmdm, 'leiden', method='wilcoxon', key_added='leiden_clust')
# result = bmdm_all.uns['rank_genes_groups']['names'].dtype.names
# cluster_diffexp = pd.DataFrame({group + '_' + key[:1]: bmdm_all.uns['rank_genes_groups'][key][group]
#                                 for group in result for key in ['names', 'pvals_adj', 'logfoldchanges']})

# writing result for each cluster to xlsx file
clust0 = sc.get.rank_genes_groups_df(bmdm, group='0', key='leiden_clust', pval_cutoff=0.05, log2fc_min=1.5)
clust1 = sc.get.rank_genes_groups_df(bmdm, group='1', key='leiden_clust', pval_cutoff=0.05, log2fc_min=1.5)
clust2 = sc.get.rank_genes_groups_df(bmdm, group='2', key='leiden_clust', pval_cutoff=0.05, log2fc_min=1.5)
clust3 = sc.get.rank_genes_groups_df(bmdm, group='3', key='leiden_clust', pval_cutoff=0.05, log2fc_min=1.5)
clust4 = sc.get.rank_genes_groups_df(bmdm, group='4', key='leiden_clust', pval_cutoff=0.05, log2fc_min=1.5)
clust5 = sc.get.rank_genes_groups_df(bmdm, group='5', key='leiden_clust', pval_cutoff=0.05, log2fc_min=1.5)

# write results to xlsx
# write results to xlsx
excelpath = '/Users/katebridges/PyCharmProjects/test/gse117176_leiden_20211105.xlsx'
writer = pd.ExcelWriter(excelpath, engine='xlsxwriter')

clust0.to_excel(writer, sheet_name='Cluster0')
clust1.to_excel(writer, sheet_name='Cluster1')
clust2.to_excel(writer, sheet_name='Cluster2')
clust3.to_excel(writer, sheet_name='Cluster3')
clust4.to_excel(writer, sheet_name='Cluster4')
clust5.to_excel(writer, sheet_name='Cluster5')

writer.save()


def color_byexpr(dat, gene, embed, embed_type, color_map, fig_size, cbar_max):
    fig, ax = plt.subplots(figsize=fig_size)
    scatter_x = embed[:, 0]
    scatter_y = embed[:, 1]
    fig1 = ax.scatter(scatter_x, scatter_y, c=np.array(dat[:, gene].X.todense()).flatten(), s=4, alpha=0.85, cmap=color_map,
                      vmin=0, vmax=cbar_max)
    plt.colorbar(fig1)
    plt.title(gene)
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])
    plt.xlabel(embed_type + '1')
    plt.ylabel(embed_type + '2')
    plt.show()


def viz_by(embed, cond, cmap_dict, embed_type, fig_size):
    fig, ax = plt.subplots(figsize=fig_size)
    scatter_x = embed[:, 0]
    scatter_y = embed[:, 1]
    b = 0
    for g in np.unique(cond):
        i = np.where(cond == g)
        ax.scatter(scatter_x[i], scatter_y[i], label=g, s=4, c=cmap_dict[b])
        b = b + 1
    ax.legend(loc='upper right', prop={'size': 8})
    plt.xlabel(embed_type + '1')
    plt.ylabel(embed_type + '2')
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])
    plt.show()


viz_by(bmdm.obsm['X_umap'], bmdm.obs['leiden'], sns.husl_palette(6), 'UMAP', (5, 5))

kmj_genes = ['Il12b', 'Il6', 'Tnf', 'Ccl5']
for j in kmj_genes:
    color_byexpr(bmdm, j, bmdm.obsm['X_umap'], 'UMAP', 'Reds', (6, 5), None)

bar_ind = np.unique(bmdm.obs['leiden'])
for gene_name in kmj_genes:
    dat = np.array(bmdm[:, gene_name].X.todense()).flatten()
    dat_stat = np.zeros((len(np.unique(bmdm.obs['leiden'])), 3))
    b = 0
    for g in np.unique(bmdm.obs['leiden']):
        i = np.where(bmdm.obs['leiden'] == g)[0]
        ci_info = bs.bootstrap(dat[i], stat_func=bs_stats.mean)
        dat_stat[b, 0] = ci_info.value
        dat_stat[b, 1] = dat_stat[b, 0] - ci_info.lower_bound
        dat_stat[b, 2] = ci_info.upper_bound - dat_stat[b, 0]
        b = b + 1
    fig, ax = plt.subplots(figsize=(4, 5))
    barlist = ax.bar(bar_ind, dat_stat[:, 0], yerr=[dat_stat[:, 1], dat_stat[:, 2]], align='center', ecolor='black', capsize=10)
    barlist[0].set_color(sns.husl_palette(6)[0])
    barlist[1].set_color(sns.husl_palette(6)[1])
    barlist[2].set_color(sns.husl_palette(6)[2])
    barlist[3].set_color(sns.husl_palette(6)[3])
    barlist[4].set_color(sns.husl_palette(6)[4])
    barlist[5].set_color(sns.husl_palette(6)[5])
    plt.title(gene_name)
    ax.set_ylabel('ln[mRNA counts + 1]')
    ax.set_xticks(np.arange(len(bar_ind)))
    ax.set_xticklabels(bar_ind)
    # ax.set_xlim([-0.5, len(bar_ind) + 0.5])
    plt.tight_layout()

# taking a look at mac/DC-specific marker expr
helft_mac = ['Ptplad2', '1810011H11Rik', 'Tlr4', 'Fgd4', 'Sqrdl', 'Csf3r', 'Plod1', 'Tom1', 'Pld3', 'Tpp1', 'Ctsd',
             'Lamp2', 'Pla2g4a', 'Fcgr1', 'Mr1', 'Mertk', 'Cd14', 'Tbxas1', 'Fcgr3', 'Sepp1', 'Cd164', 'Tcn2', 'Dok3',
             'Ctsl', 'Tspan14', 'Itgam', 'Adgre1', 'Sirpa', 'Csf1r']
# most of above genes were nonvariable (?) and aren't included in this version of the data
helft_dc = ['Adam19', 'Ccr7', 'Gpr132', 'H2-Eb2', 'Hmgn3', 'Kit', 'Klri1', 'Kmo', 'P2ry10', 'Pvrl1', 'Rab30',
            'Slamf7', 'Traf1', 'Zbtb46', 'Dpp4', 'Runx3', 'Itgax', 'Cd24a', 'Flt3', 'Ly75', 'Pdcd1lg2', 'Irf4',
            'H2-Ab1', 'H2-Eb1', 'Cd74', 'H2-Ob', 'H2-Aa', 'Cd86', 'Cd80', 'H2-DMb2']

mac_upd = []
dc_upd = []

for i in helft_mac:
    g = np.where(bmdm.var_names == i)[0]
    if g.size:
        mac_upd.append(i)

for f in dc_upd:
    color_byexpr(bmdm, f, bmdm.obsm['X_umap'], 'UMAP', 'Reds', (6, 5), None)

# scoring by entire signature
sc.tl.score_genes(bmdm, mac_upd, score_name='UCG_Mac_score')
sc.tl.score_genes(bmdm, dc_upd, score_name='UCG_DC_score')


def color_byscore(score, embed, embed_type, color_map, fig_size, title):
    fig, ax = plt.subplots(figsize=fig_size)
    scatter_x = embed[:, 0]
    scatter_y = embed[:, 1]
    fig1 = ax.scatter(scatter_x, scatter_y, c=score, s=4, cmap=color_map)
    plt.colorbar(fig1)
    plt.title(title)
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])
    plt.xlabel(embed_type + '1')
    plt.ylabel(embed_type + '2')
    plt.show()


color_byscore(bmdm.obs['UCG_Mac_score'], bmdm.obsm['X_umap'], 'UMAP', 'RdYlBu_r', (6, 5), 'UCG Macrophage score')

UCG_all = mac_upd + dc_upd
all_ucg = np.zeros((len(UCG_all), len(np.unique(bmdm.obs['leiden']))))

for j in np.unique(bmdm.obs['leiden']):
    dat = bmdm[np.where(bmdm.obs['leiden'] == j)[0], :]
    for k in range(len(UCG_all)):
        all_ucg[k, int(j)] = np.mean(dat[:, UCG_all[k]].X)

ucg_labels = np.concatenate((np.tile(0, len(mac_upd)), np.tile(1, len(dc_upd))))

network_pal = sns.husl_palette(len(np.unique(ucg_labels)))
network_lut = dict(zip(np.unique(ucg_labels), network_pal))
network_colors = pd.Series(ucg_labels).map(network_lut)

cg = sns.clustermap(all_ucg, yticklabels=UCG_all, xticklabels=np.arange(6), cmap='RdYlBu_r', linewidths=0.1,
               linecolor='black', rasterized=False, col_cluster=True, row_cluster=True, center=0, z_score=0,
               row_colors=network_colors.values, figsize=(6, 12))
cg.savefig('gse_heatmap.png')


# only really seeing "DC-like" cluster in M1 condition so perhaps we limit analysis?
bmdm_m1 = bmdm[np.where(bmdm.obs['sample'] == 'M1')[0], :]
del bmdm_m1.obsm, bmdm_m1.obsp, bmdm_m1.uns

sc.tl.pca(bmdm_m1, svd_solver='auto')
sc.pl.pca(bmdm_m1, color='sample')
sc.pp.neighbors(bmdm_m1)  # using with default parameters

# UMAP visualization
sc.tl.umap(bmdm_m1)  # default initial position
sc.tl.leiden(bmdm_m1, resolution=0.2)

sc.tl.rank_genes_groups(bmdm_m1, 'leiden', method='wilcoxon', key_added='leiden_clust')
# result = bmdm_all.uns['rank_genes_groups']['names'].dtype.names
# cluster_diffexp = pd.DataFrame({group + '_' + key[:1]: bmdm_all.uns['rank_genes_groups'][key][group]
#                                 for group in result for key in ['names', 'pvals_adj', 'logfoldchanges']})

# writing result for each cluster to xlsx file
clust0 = sc.get.rank_genes_groups_df(bmdm_m1, group='0', key='leiden_clust', pval_cutoff=0.05, log2fc_min=1.5)
clust1 = sc.get.rank_genes_groups_df(bmdm_m1, group='1', key='leiden_clust', pval_cutoff=0.05, log2fc_min=1.5)
clust2 = sc.get.rank_genes_groups_df(bmdm_m1, group='2', key='leiden_clust', pval_cutoff=0.05, log2fc_min=1.5)

# write results to xlsx
# write results to xlsx
excelpath = '/Users/katebridges/PyCharmProjects/test/gse117176_leiden_20211105_M1ONLY.xlsx'
writer = pd.ExcelWriter(excelpath, engine='xlsxwriter')

clust0.to_excel(writer, sheet_name='Cluster0')
clust1.to_excel(writer, sheet_name='Cluster1')
clust2.to_excel(writer, sheet_name='Cluster2')

writer.save()

viz_by(bmdm_m1.obsm['X_umap'], bmdm_m1.obs['leiden'], sns.husl_palette(3), 'UMAP', (5, 5))
color_byexpr(bmdm_m1, 'Ccl22', bmdm_m1.obsm['X_umap'], 'UMAP', 'Reds', (6, 5), None)

for f in dc_upd:
    color_byexpr(bmdm_m1, f, bmdm_m1.obsm['X_umap'], 'UMAP', 'Reds', (6, 5), None)

color_byexpr(bmdm_m1, 'Tcn2', bmdm_m1.obsm['X_umap'], 'UMAP', 'Reds', (6, 5), None)

# scoring by entire signature
sc.tl.score_genes(bmdm_m1, mac_upd, score_name='UCG_Mac_score')
sc.tl.score_genes(bmdm_m1, dc_upd, score_name='UCG_DC_score')

color_byscore(bmdm_m1.obs['UCG_DC_score'], bmdm_m1.obsm['X_umap'], 'UMAP', 'RdYlBu_r', (6, 5), 'UCG DC score (M1)')

m1_ucg = np.zeros((len(UCG_all), len(np.unique(bmdm_m1.obs['leiden']))))

for j in np.unique(bmdm_m1.obs['leiden']):
    dat = bmdm_m1[np.where(bmdm_m1.obs['leiden'] == j)[0], :]
    for k in range(len(UCG_all)):
        m1_ucg[k, int(j)] = np.mean(dat[:, UCG_all[k]].X)

cg = sns.clustermap(m1_ucg, yticklabels=UCG_all, xticklabels=np.arange(3), cmap='RdYlBu_r', linewidths=0.1,
               linecolor='black', rasterized=False, col_cluster=True, row_cluster=True, center=0, z_score=0,
               row_colors=network_colors.values, figsize=(6, 12))
cg.savefig('gse_heatmap_m1.png')
