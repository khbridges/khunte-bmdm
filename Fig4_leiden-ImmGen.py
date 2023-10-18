import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import scvelo as scv
import loompy

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

# reading from backup
bmdm_fig3 = sc.read('/Users/katebridges/Downloads/bmdm_object_20220509.h5ad')

loom_dir = '/Users/katebridges/Downloads/AMR-LOOM/'
m1_ldata = scv.read(loom_dir + 'M1-H_MMT.loom', cache=False)
loompy.combine([loom_dir + 'M0-H_MMT.loom', loom_dir + 'M1-H_MMT.loom', loom_dir + 'M2-H_MMT.loom'],
               output_file=loom_dir + 'merged-BMDM.loom')
merged_ldata = scv.read(loom_dir + 'merged-BMDM.loom', cache=False)

# integrating with loom files for RNA velo
bmdm_fig3 = scv.utils.merge(bmdm_fig3, merged_ldata)

# running basic RNA velocity analyses
scv.pl.proportions(bmdm_fig3, groupby='sample')
scv.tl.velocity(bmdm_fig3)
scv.tl.velocity_graph(bmdm_fig3)
scv.pl.velocity_embedding_stream(bmdm_fig3, basis='umap', color='leiden')

# PAGA for reiteration of cross-condition cluster similarity for clust 5?
bmdm_fig3.obs['sample_leiden'] = [bmdm_fig3.obs['sample'][h] + ' ' + bmdm_fig3.obs['leiden'][h] for h in bmdm_fig3.obs_names]
# bmdm_fig3_plotting obj - remove 'sample_leiden' categories with 1 cell
sc.tl.paga(bmdm_fig3, groups='sample_leiden', use_rna_velocity=False)

paga_palette = {'M0 0': sns.husl_palette(6)[0], 'M0 1': sns.husl_palette(6)[1], 'M0 4': sns.husl_palette(6)[4], 'M0 5': sns.husl_palette(6)[5], 'M1 2': sns.husl_palette(6)[2], 'M1 3': sns.husl_palette(6)[3], 'M1 5': sns.husl_palette(6)[5], 'M2 0': sns.husl_palette(6)[0], 'M2 1': sns.husl_palette(6)[1],  'M2 4': sns.husl_palette(6)[4], 'M2 5': sns.husl_palette(6)[5]}
bmdm_fig3_plotting.uns['sample_leiden_colors'] = list(paga_palette.values())
sc.tl.paga(bmdm_fig3_plotting, groups='sample_leiden', use_rna_velocity=False)
sc.pl.paga(bmdm_fig3_plotting, threshold=0.1, node_size_scale=5, node_size_power=0.4)

# get DEGs for cluster of all data
sc.tl.rank_genes_groups(bmdm_fig3, 'leiden', method='wilcoxon', key_added='leiden_clust')

# writing result for each cluster to xlsx file
clust0 = sc.get.rank_genes_groups_df(bmdm_fig3, group='0', key='leiden_clust', pval_cutoff=0.05, log2fc_min=1.5)
clust1 = sc.get.rank_genes_groups_df(bmdm_fig3, group='1', key='leiden_clust', pval_cutoff=0.05, log2fc_min=1.5)
clust2 = sc.get.rank_genes_groups_df(bmdm_fig3, group='2', key='leiden_clust', pval_cutoff=0.05, log2fc_min=1.5)
clust3 = sc.get.rank_genes_groups_df(bmdm_fig3, group='3', key='leiden_clust', pval_cutoff=0.05, log2fc_min=1.5)
clust4 = sc.get.rank_genes_groups_df(bmdm_fig3, group='4', key='leiden_clust', pval_cutoff=0.05, log2fc_min=1.5)
clust5 = sc.get.rank_genes_groups_df(bmdm_fig3, group='5', key='leiden_clust', pval_cutoff=0.05, log2fc_min=1.5)

# write results to xlsx
excelpath = '/Users/katebridges/PyCharmProjects/test/bmdm_leiden_20211019.xlsx'
writer = pd.ExcelWriter(excelpath, engine='xlsxwriter')

clust0.to_excel(writer, sheet_name='Cluster0')
clust1.to_excel(writer, sheet_name='Cluster1')
clust2.to_excel(writer, sheet_name='Cluster2')
clust3.to_excel(writer, sheet_name='Cluster3')
clust4.to_excel(writer, sheet_name='Cluster4')
clust5.to_excel(writer, sheet_name='Cluster5')

writer.save()

# UMAP VIZ by leiden cluster
leiden_dict = {np.unique(bmdm_fig3.obs['leiden'])[m]: sns.husl_palette(6)[m] for m in np.arange(len(np.unique(bmdm_fig3.obs['leiden'])))}
sc.pl.umap(bmdm_fig3, color='leiden', palette=leiden_dict)

# HEATMAP for leiden clusters x mac/DC signature genes (sig genes from Helft et al.)
helft_mac = ['Ptplad2', '1810011H11Rik', 'Tlr4', 'Fgd4', 'Sqrdl', 'Csf3r', 'Plod1', 'Tom1', 'Pld3', 'Tpp1', 'Ctsd',
             'Lamp2', 'Pla2g4a', 'Fcgr1', 'Mr1', 'Mertk', 'Cd14', 'Tbxas1', 'Fcgr3', 'Sepp1', 'Cd164', 'Tcn2', 'Dok3',
             'Ctsl', 'Tspan14', 'Itgam', 'Adgre1', 'Sirpa', 'Csf1r']

helft_dc = ['Adam19', 'Ccr7', 'Gpr132', 'H2-Eb2', 'Hmgn3', 'Kit', 'Klri1', 'Kmo', 'P2ry10', 'Pvrl1', 'Rab30',
            'Slamf7', 'Traf1', 'Zbtb46', 'Dpp4', 'Runx3', 'Itgax', 'Cd24a', 'Flt3', 'Ly75', 'Pdcd1lg2', 'Irf4',
            'H2-Ab1', 'H2-Eb1', 'Cd74', 'H2-Ob', 'H2-Aa', 'Cd86', 'Cd80', 'H2-DMb2']
mac_upd = []
dc_upd = []

# to sort out genes that were nonvariable/excluded during preprocessing:
for j in helft_mac:
    h = np.where(bmdm_all.var_names == j)[0]
    if h.size:
        mac_upd.append(j)
for i in helft_dc:
    g = np.where(bmdm_all.var_names == i)[0]
    if g.size:
        dc_upd.append(i)

UCG_all = mac_upd + dc_upd
all_ucg = np.zeros((len(UCG_all), len(np.unique(bmdm_fig3.obs['leiden']))))

for j in np.unique(bmdm_all.obs['leiden']):
    dat = bmdm_fig3[np.where(bmdm_fig3.obs['leiden'] == j)[0], :]
    for k in range(len(UCG_all)):
        all_ucg[k, int(j)] = np.mean(dat[:, UCG_all[k]].X)

# plotting avg expr of signature genes for each leiden cluster
# row colors to designate whether signature belongs to mac or DC
ucg_labels = np.concatenate((np.tile(0, len(mac_upd)), np.tile(1, len(dc_upd))))

network_pal = sns.husl_palette(len(np.unique(ucg_labels)))
network_lut = dict(zip(np.unique(ucg_labels), network_pal))
network_colors = pd.Series(ucg_labels).map(network_lut)

cg = sns.clustermap(all_ucg, yticklabels=UCG_all, xticklabels=np.arange(6), cmap='RdYlBu_r', linewidths=0.1,
                    linecolor='black', rasterized=False, col_cluster=True, row_cluster=True, center=0, z_score=0,
                    row_colors=network_colors.values, figsize=(6, 12))
cg.savefig('fig3_heatmap.png')

# IMMGEN ANALYSES
# reading in ImmGen samples
mnp = pd.read_csv('/Users/katebridges/Downloads/GSE122108_Gene_count_table.csv')

# finding unique samples (helps to identify replicates for combining)
dc_samples = mnp.columns[mnp.columns.str.contains('DC')]
# dc_unique = []
# for h in dc_samples:
#     dc_unique.append(h[:-2])

# reading in sample code conversions (to simplified groupings)
samplecode = pd.read_csv('/Users/katebridges/PycharmProjects/khunte-bmdm/ImmGen_DC_samplecodes.csv')

org_lab = [None] * len(dc_samples)
b = 0
for g in samplecode['sample_code'].values:
    i = np.where(dc_samples.str.contains(g))[0]
    for h in i:
        org_lab[h] = samplecode['simplified'][b]
    b = b + 1

# restricting mnp dataset to 12k genes included in preprocessed scRNA-seq - need to read in dataset for future ref
genes = np.intersect1d(mnp['GeneSymbol'].values, bmdm_fig3.var_names.values)
gene_ind = np.zeros(len(genes))
k = 0
for j in genes:
    gene_ind[k] = np.where(mnp['GeneSymbol'] == j)[0][0]
    k = k + 1
gene_ind = gene_ind.astype('int')
mnp_restricted = mnp.iloc[gene_ind, :]

# finding averages gene expr counts
dc_simplified = np.zeros((np.unique(samplecode['simplified']).shape[0], len(genes)))
w = 0
for q in np.unique(samplecode['simplified']):
    i = np.where(samplecode['simplified'] == q)[0]
    if len(i) > 1:
        g = np.array([])
        for m in i:
            g = np.concatenate((g, mnp.columns[mnp.columns.str.contains(samplecode['sample_code'][m])]))
        dc_simplified[w] = (mnp_restricted[g].sum(axis=1))/(1000*len(g))
    else:
        dc_simplified[w] = (mnp_restricted[mnp.columns[mnp.columns.str.contains(samplecode['sample_code'][i[0]])]].sum(axis=1))/(1000*len(mnp.columns[mnp.columns.str.contains(samplecode['sample_code'][i[0]])]))
    w = w + 1

# now let's create a pandas DataFrame with this data
dc_immgen = pd.DataFrame(data=dc_simplified, index=np.unique(samplecode['simplified']), columns=genes)

# now to create a similar dataframe across 6 clusters in bmdm data
bmdm_counts = bmdm_fig3[:, genes].layers['counts'].todense()
bmdm_simplified = np.zeros((len(np.unique(bmdm_fig3.obs['leiden'])), len(genes)))
for d in np.unique(bmdm_fig3.obs['leiden']):
    i = np.where(bmdm_fig3.obs['leiden'] == d)[0]
    bmdm_simplified[int(d), :] = np.mean(bmdm_counts[i, :], axis=0)

cluster_name = np.array(['0 (Ctrl)', '1 (+IL-4)', '2 (+LPS+IFNg)', '3 (+LPS+IFNg)', '4 (Ctrl)', '5 (DC-like)'])
dc_bmdm = pd.DataFrame(data=bmdm_simplified, index=cluster_name, columns=genes)

# now to combine dataframes, limit to most UCG signature genes
dc = pd.concat([dc_bmdm, dc_immgen])
top_ind = dc.columns[np.argsort(dc.std())][(dc.shape[1] - int(dc.shape[1]*0.5)):]
ig_ind = np.intersect1d(UCG_all, genes)

ucg_labels = np.array([1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1,
                       1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1])

network_pal = sns.husl_palette(len(np.unique(ucg_labels)))
network_lut = dict(zip(np.unique(ucg_labels), network_pal))
network_colors = pd.Series(ucg_labels).map(network_lut)

# visualize
ig = sns.clustermap(np.log(dc[ig_ind]+1), yticklabels=dc.index, cmap='viridis_r', linecolor='black', rasterized=False,
                    col_cluster=True, row_cluster=True, xticklabels=ig_ind, col_colors=network_colors.values)
ig.savefig('immgen_dc.png')


# ENRICHMENT for DC subsets
dc_geneset = {'mregDC': ['Cd80', 'Cd40', 'Cd83', 'Cd86', 'Relb', 'Cd274', 'Pdcd1lg2', 'Cd200', 'Fas', 'Socs1',
                         'Socs2', 'Aldh1a2', 'Ccr7', 'Fscn1', 'Il4ra', 'Il4i1', 'Myo1g', 'Cxcl16', 'Adam8', 'Icam1',
                         'Marcks', 'Marcksl1'],
              'DC1': ['Xcr1', 'Clec9a', 'Cadm1', 'Naaa'],
              'DC2': ['Itgam', 'Cd209a', 'Sirpa', 'H2-DMb2']}

sc.tl.score_genes(bmdm_fig3, dc_geneset['mregDC'], score_name='mregDC_score')
sc.tl.score_genes(bmdm_fig3, dc_geneset['DC1'], score_name='DC1_score')
sc.tl.score_genes(bmdm_fig3, dc_geneset['DC2'], score_name='DC2_score')
