#!/usr/bin/env python
# coding: utf-8


import sys
import os
import pandas as pd
import pwseqdist as pw
import tcrdist
from tcrdist.repertoire import TCRrep
from hybrid.feasible import tcr_ab_to_network
from hybrid.feasible import get_hla_dictionaries
from hybrid.feasible import compute_hla_feasibility
from tcrdist.public import _neighbors_fixed_radius
import pandas as pd
import networkx as nx
import community.community_louvain as community_louvain
from tcrdist.repertoire import TCRrep
import networkx as nx
from matplotlib import pyplot as plt
import random
import matplotlib.pyplot as plt
import networkx as nx
import random


# # Data

# In[41]:


data_path =  './data/'
# Load a compiled file containing all phenotyped cells
clones_file = pd.read_csv(os.path.join(data_path, '1_Tcell_tcr_meta_filter.csv'))
clones_file.columns


# In[42]:


clones_file.drop(columns=['Unnamed: 0'], inplace=True)
print(clones_file['ct_level_4'].unique())
print(len(clones_file))


# # TCRrep and sequence network

# In[36]:


def Tcrrep_and_sequance_network(clones_file, targetcell,edge_threshold,filename,resolution=1):
    #df = clones_file[(clones_file['ct_level_4'] == targetcell)]
    df = clones_file[clones_file['ct_level_4'].isin(targetcell)]
    df['v_b_gene'] = df['v_b_gene'].apply(lambda x : f"{x}*01")
    df['j_b_gene'] = df['j_b_gene'].apply(lambda x : f"{x}*01")
    clone_cols = ['Condition', 'cdr3_b_aa',
              'v_b_gene','j_b_gene','ct_level_4','barcode','Donor']
    tr = TCRrep(cell_df = df[clone_cols], chains = ['beta'], organism = 'human')
     # 数量不能超过1w
   
    network = list()
    for i,n in enumerate(_neighbors_fixed_radius(tr.pw_beta, edge_threshold)):
      for j in n:
          if i != j:
              network.append((
                  i,                                 # 'node_1' - row index
                  j,                                 # 'node_2' - column index
                  (tr.pw_beta)[i,j],          # 'dist'- gets the distance between TCR(i,j)
                  tr.clone_df['v_b_gene'].iloc[i],   # 'v_b_gene_1' - v beta gene of clone i
                  tr.clone_df['v_b_gene'].iloc[j],   # 'v_b_gene_2' - v beta gene of clone j
                  tr.clone_df['cdr3_b_aa'].iloc[i],  # 'cdr3_b_aa_1' - cdr3 beta of clone i
                  tr.clone_df['cdr3_b_aa'].iloc[j],  # 'cdr3_b_aa_2' - cdr3 beta of clone j
                  tr.clone_df['count'].iloc[i],      #  count 
                  tr.clone_df['count'].iloc[j],      #  count
                  tr.clone_df['Condition'].iloc[i],      #  Condition 
                  tr.clone_df['Condition'].iloc[j],      #  Condition
                  len(n)-1))                         # 'K_neighbors' - number of neighbors
    
    cols = ['node_1', 'node_2', 'dist', 
          'v_b_gene_1', 'v_b_gene_2', 
          'cdr3_b_aa_1','cdr3_b_aa_2', 
          'count_1', 'count_2',
          'Condition_1', 'Condition_2',
          'K_neighbors']
    # Store the <network> edge list as a DataFrame.        
    df_net = pd.DataFrame(network, columns = cols)
    df_net['weight'] = (edge_threshold - df_net['dist'])/edge_threshold
    G = nx.from_pandas_edgelist(pd.DataFrame({'source' : df_net['node_1'],'target' : df_net['node_2'], 'weight' :df_net['weight']}))
    partition= community_louvain.best_partition(G,resolution=resolution) # resolution默认为1,较高的值能导致更小、数量更多的社区
    # Change partition such that cluster Id is in descending order based on community size 
    partitions_by_cluster_size = list(pd.Series(partition.values()).value_counts().index)
    partition_reorder = {id:rank for id,rank in zip(partitions_by_cluster_size, 
                                                    range(len(partitions_by_cluster_size)))}
    partition = {k:partition_reorder.get(v) for k,v in partition.items()}
    
    from tcrdist.html_colors import get_html_colors
    clusters = [i for i in pd.Series(partition.values()).value_counts().index]
    colors   = get_html_colors(len(clusters))
    cluster_to_color = {cluster:color for cluster,color, in zip(clusters,colors)}
    print(G)
    print(nx.number_connected_components(G), "connected components")
    plt.figure(1, figsize=(16, 16))
    # layout graphs with positions using graphviz neato
    pos = nx.nx_agraph.graphviz_layout(G, prog="neato") # neato
    
    # color nodes the same in each connected subgraph
    C = (G.subgraph(c) for c in nx.connected_components(G))
    ns = df_net['node_1'].to_list() + df_net['node_2'].to_list()
    cts = df_net['count_1'].to_list() + df_net['count_2'].to_list()
    node_cts = {node: x for node, x in zip(ns,cts)}
    # for g in C:
    #     out_dict =dict()
    #     count = 0
    #     v_so_far = dict()
    #     for k,v in {i : partition.get(i) for i in g.nodes}.items():
    #         if v <= 100 and v_so_far.get(v, True):
    #             out_dict[k] =v 
    #         else:
    #             out_dict[k] = ""
    #         v_so_far[v] = False
    #    #c = [random.random()] * nx.number_of_nodes(g)  # random color...
    #     nx.draw(g, pos, alpha= .5, node_size = [5*node_cts.get(x) for x in g.nodes], 
    #             node_color = [cluster_to_color.get(partition.get(i)) for i in g.nodes],
    #             edge_color = "gray",
    #             labels = out_dict,
    #             vmin=0.0, vmax=1.0, with_labels=True, font_size = 20, font_color = "violet")

    # 可视化
    for g in C:
        out_dict = dict()
        count = 0
        v_so_far = dict()
        for k, v in {i: partition.get(i) for i in g.nodes}.items():
            if v <= 100 and v_so_far.get(v, True):
                out_dict[k] = v 
            else:
                out_dict[k] = ""
            v_so_far[v] = False
    
        nx.draw(g, pos, alpha=.5, node_size=[5 * node_cts.get(x) for x in g.nodes], 
                node_color=[cluster_to_color.get(partition.get(i)) for i in g.nodes],
                edge_color="gray",
                labels=out_dict,
                vmin=0.0, vmax=1.0, with_labels=True, font_size=30, font_color="black", font_weight='bold')
    
    plt.savefig(f"figures/1_{filename}_edge_{edge_threshold}.pdf")
    
    tr.clone_df['cluster_id'] = [str(partition.get(i)) if partition.get(i) is not None else None for i in tr.clone_df.index]
    tr.clone_df['cluster_id_cat'] =  pd.Series(tr.clone_df['cluster_id'].to_list(), dtype="category").cat.reorder_categories(tr.clone_df['cluster_id'].value_counts().index)
    
    tr.clone_df.to_csv(f"data/2_{filename}_edge_{edge_threshold}.csv")
    return tr.clone_df


# In[47]:


def Tcrrep_and_sequance_network(clones_file, targetcell,edge_threshold,filename,resolution=1):
    #df = clones_file[(clones_file['ct_level_4'] == targetcell)]
    df = clones_file[clones_file['ct_level_4'].isin(targetcell)]
    df['v_b_gene'] = df['v_b_gene'].apply(lambda x : f"{x}*01")
    df['j_b_gene'] = df['j_b_gene'].apply(lambda x : f"{x}*01")
    clone_cols = ['Condition', 'cdr3_b_aa',
              'v_b_gene','j_b_gene','ct_level_4','barcode','Donor']
    tr = TCRrep(cell_df = df[clone_cols], chains = ['beta'], organism = 'human')
     # 数量不能超过1w
   
    network = list()
    for i,n in enumerate(_neighbors_fixed_radius(tr.pw_beta, edge_threshold)):
      for j in n:
          if i != j:
              network.append((
                  i,                                 # 'node_1' - row index
                  j,                                 # 'node_2' - column index
                  (tr.pw_beta)[i,j],          # 'dist'- gets the distance between TCR(i,j)
                  tr.clone_df['v_b_gene'].iloc[i],   # 'v_b_gene_1' - v beta gene of clone i
                  tr.clone_df['v_b_gene'].iloc[j],   # 'v_b_gene_2' - v beta gene of clone j
                  tr.clone_df['cdr3_b_aa'].iloc[i],  # 'cdr3_b_aa_1' - cdr3 beta of clone i
                  tr.clone_df['cdr3_b_aa'].iloc[j],  # 'cdr3_b_aa_2' - cdr3 beta of clone j
                  tr.clone_df['count'].iloc[i],      #  count 
                  tr.clone_df['count'].iloc[j],      #  count
                  tr.clone_df['Condition'].iloc[i],      #  Condition 
                  tr.clone_df['Condition'].iloc[j],      #  Condition
                  len(n)-1))                         # 'K_neighbors' - number of neighbors
    
    cols = ['node_1', 'node_2', 'dist', 
          'v_b_gene_1', 'v_b_gene_2', 
          'cdr3_b_aa_1','cdr3_b_aa_2', 
          'count_1', 'count_2',
          'Condition_1', 'Condition_2',
          'K_neighbors']
    # Store the <network> edge list as a DataFrame.        
    df_net = pd.DataFrame(network, columns = cols)
    df_net['weight'] = (edge_threshold - df_net['dist'])/edge_threshold
    G = nx.from_pandas_edgelist(pd.DataFrame({'source' : df_net['node_1'],'target' : df_net['node_2'], 'weight' :df_net['weight']}))
    partition= community_louvain.best_partition(G,resolution=resolution) # resolution默认为1,较高的值能导致更小、数量更多的社区
    # Change partition such that cluster Id is in descending order based on community size 
    partitions_by_cluster_size = list(pd.Series(partition.values()).value_counts().index)
    partition_reorder = {id:rank for id,rank in zip(partitions_by_cluster_size, 
                                                    range(len(partitions_by_cluster_size)))}
    partition = {k:partition_reorder.get(v) for k,v in partition.items()}
    
    from tcrdist.html_colors import get_html_colors
    clusters = [i for i in pd.Series(partition.values()).value_counts().index]
    colors   = get_html_colors(len(clusters))
    cluster_to_color = {cluster:color for cluster,color, in zip(clusters,colors)}
    print(G)
    print(nx.number_connected_components(G), "connected components")
    plt.figure(1, figsize=(16, 16))
    # layout graphs with positions using graphviz neato
    pos = nx.nx_agraph.graphviz_layout(G, prog="neato") # neato
    
    # color nodes the same in each connected subgraph
    C = (G.subgraph(c) for c in nx.connected_components(G))
    ns = df_net['node_1'].to_list() + df_net['node_2'].to_list()
    cts = df_net['count_1'].to_list() + df_net['count_2'].to_list()
    node_cts = {node: x for node, x in zip(ns,cts)}
    for g in C:
        out_dict =dict()
        count = 0
        v_so_far = dict()
        for k,v in {i : partition.get(i) for i in g.nodes}.items():
            if v <= 100 and v_so_far.get(v, True):
                out_dict[k] =v 
            else:
                out_dict[k] = ""
            v_so_far[v] = False
       #c = [random.random()] * nx.number_of_nodes(g)  # random color...
        nx.draw(g, pos, alpha= .5, node_size = [5*node_cts.get(x) for x in g.nodes], 
                node_color = [cluster_to_color.get(partition.get(i)) for i in g.nodes],
                edge_color = "gray",
                labels = out_dict,
                vmin=0.0, vmax=1.0, with_labels=True, font_size = 20, font_color = "violet")

 
    
    plt.savefig(f"figures/1_{filename}_edge_{edge_threshold}.pdf")
    
    tr.clone_df['cluster_id'] = [str(partition.get(i)) if partition.get(i) is not None else None for i in tr.clone_df.index]
    tr.clone_df['cluster_id_cat'] =  pd.Series(tr.clone_df['cluster_id'].to_list(), dtype="category").cat.reorder_categories(tr.clone_df['cluster_id'].value_counts().index)
    
    tr.clone_df.to_csv(f"data/2_{filename}_edge_{edge_threshold}.csv")
    return tr.clone_df


# In[48]:


clones_file['ct_level_4'].value_counts()


# In[49]:


#已改源代码把10000细胞限制改成40000


# In[51]:


CD8T_all_res =Tcrrep_and_sequance_network( clones_file = clones_file,
                            targetcell = ["CD8_Tem","CD8_Temra",
                                          "CD8_NELL2","CD8_Tn"],
                            edge_threshold = 35,
                           filename = 'CD8T_all' )


# In[ ]:


CD8T_all_res =Tcrrep_and_sequance_network( clones_file = clones_file,
                            targetcell = ["CD8_Tem","CD8_Temra",
                                          "CD8_NELL2","CD8_Tn"],
                            edge_threshold = 25,
                           filename = 'CD8T_all' )


# In[24]:


CD8T_all_res =Tcrrep_and_sequance_network( clones_file = clones_file,
                            targetcell = ["CD8_Tem","CD8_Temra",
                                          "CD8_NELL2","CD8_Tn"],
                            edge_threshold = 15,
                           filename = 'CD8T_all' )


# In[25]:


CD8T_all_res =Tcrrep_and_sequance_network( clones_file = clones_file,
                            targetcell = ["CD8_Tem","CD8_Temra",
                                          "CD8_NELL2","CD8_Tn"],
                            edge_threshold = 50,
                           filename = 'CD8T_all' )


# In[26]:


CD8T_all_res =Tcrrep_and_sequance_network( clones_file = clones_file,
                            targetcell = ["CD8_Tem","CD8_Temra",
                                          "CD8_NELL2","CD8_Tn"],
                            edge_threshold = 60,
                           filename = 'CD8T_all' )


# In[27]:


CD8T_all_res =Tcrrep_and_sequance_network( clones_file = clones_file,
                            targetcell = ["CD8_Tem","CD8_Temra",
                                          "CD8_NELL2","CD8_Tn"],
                            edge_threshold = 80,
                           filename = 'CD8T_all' )


# In[52]:


CD4T_all_res = Tcrrep_and_sequance_network( clones_file = clones_file,
                            targetcell = ["CD4_Tem","CD4_Temra",
                                          "CD4_Th2","CD4_T_CCR6",
                                          "CD4_Tcm","CD4_Tn","CD4_Treg"],
                            edge_threshold = 35,
                           filename = 'CD4T_all' )


# In[28]:


CD4T_all_res = Tcrrep_and_sequance_network( clones_file = clones_file,
                            targetcell = ["CD4_Tem","CD4_Temra",
                                          "CD4_Th2","CD4_T_CCR6",
                                          "CD4_Tcm","CD4_Tn","CD4_Treg"],
                            edge_threshold = 25,
                           filename = 'CD4T_all' )


# In[29]:


CD4T_all_res = Tcrrep_and_sequance_network( clones_file = clones_file,
                            targetcell = ["CD4_Tem","CD4_Temra",
                                          "CD4_Th2","CD4_T_CCR6",
                                          "CD4_Tcm","CD4_Tn","CD4_Treg"],
                            edge_threshold = 15,
                           filename = 'CD4T_all' )


# In[30]:


CD4T_all_res = Tcrrep_and_sequance_network( clones_file = clones_file,
                            targetcell = ["CD4_Tem","CD4_Temra",
                                          "CD4_Th2","CD4_T_CCR6",
                                          "CD4_Tcm","CD4_Tn","CD4_Treg"],
                            edge_threshold = 50,
                           filename = 'CD4T_all' )


# In[31]:


CD4T_all_res = Tcrrep_and_sequance_network( clones_file = clones_file,
                            targetcell = ["CD4_Tem","CD4_Temra",
                                          "CD4_Th2","CD4_T_CCR6",
                                          "CD4_Tcm","CD4_Tn","CD4_Treg"],
                            edge_threshold = 60,
                           filename = 'CD4T_all' )


# In[ ]:


CD4T_all_res = Tcrrep_and_sequance_network( clones_file = clones_file,
                            targetcell = ["CD4_Tem","CD4_Temra",
                                          "CD4_Th2","CD4_T_CCR6",
                                          "CD4_Tcm","CD4_Tn","CD4_Treg"],
                            edge_threshold = 80,
                           filename = 'CD4T_all' )




