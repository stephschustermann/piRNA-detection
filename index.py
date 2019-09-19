import pandas as pd

readFeaturesDataFrame1 = pd.read_csv('./featuresData/allReadFeaturesData.txt')
# use the the sequences with the features in ann
X = readFeaturesDataFrame1.iloc[:, 1:9].values

from sklearn.preprocessing import StandardScaler
X = StandardScaler().fit_transform(X)

# import HDBSCAN
import hdbscan
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
plt.switch_backend('agg')

clusterer = hdbscan.HDBSCAN(min_cluster_size=100).fit(X)

#condesed tree plot with circles
fig1, ax1 = plt.subplots()
ax1 = clusterer.condensed_tree_.plot(select_clusters=True,
                               selection_palette=sns.color_palette('deep', 8))
ax1.set
fig1.savefig('condensedTree1.png', tight_layout=True)


#condesed tree plot with out circles
fig2, ax2 = plt.subplots()
ax2 = clusterer.condensed_tree_.plot()
ax2.set
fig2.savefig('condensedTree2.png', tight_layout=True)

# heatmap the original table
fig3, ax3 = plt.subplots()
ax3 = sns.heatmap(X)
ax3.set
fig3.savefig('allDataHeatMap.png', tight_layout=True)

clusterPandas = clusterer.condensed_tree_.to_pandas()
clusterPandas.to_csv('./featuresData/condensed_tree_.txt', index=False)

# heatmap the clusterer
fig4, ax4 = plt.subplots()
ax4 = sns.heatmap(clusterPandas)
ax4.set
fig4.savefig('clustererDataHeatMap.png', tight_layout=True)

# single_linkage_tree_
fig5, ax5 = plt.subplots()
ax5 = clusterer.single_linkage_tree_.plot(cmap='viridis', colorbar=True)
ax5.set
fig5.savefig('single_linkage_tree_.png', tight_layout=True)


# minimum_spanning_tree_ let be the last one because it may trigger an exception
# due to optimized algorithm variation that skip explicit generation of the
# the spaning tree
fig6, ax6 = plt.subplots()
ax6 = clusterer.minimum_spanning_tree_.plot(edge_cmap='viridis',
                                      edge_alpha=0.6,
                                      node_size=80,
                                      edge_linewidth=2)
ax6.set
fig6.savefig('minimum_spanning_tree_.png', tight_layout=True)
