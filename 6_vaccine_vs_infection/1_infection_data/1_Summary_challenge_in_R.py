#!/usr/bin/env python
# coding: utf-8

# In[3]:


import anndata as ad
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import scanpy as sc
import numpy as np
import scanpy as sc
import pandas as pd
#from scipy.io import mmwrite


# In[18]:


adata = ad.read_h5ad("data/challenge_pbmc_cellxgene_230223.h5ad")
print(adata)


# In[27]:


# 检查数据是否存在异常
print(adata.X.shape)  # 验证矩阵大小
print(adata.obs.shape, adata.var.shape)  # 验证元数据和基因特征
print(adata.obs_names[:10])  # 检查细胞名称
print(adata.var_names[:10])  # 检查基因名称


# In[26]:


# 检查行名和索引是否一致
if adata.obs_names.equals(adata.obs.index):
    print("Row names and metadata match.")
else:
    print("Mismatch detected.")

# 检查列名和索引是否一致
if adata.var_names.equals(adata.var.index):
    print("Column names and feature metadata match.")
else:
    print("Mismatch detected.")

# 检查行名和列名唯一性
print("Unique row names:", adata.obs_names.is_unique)
print("Unique column names:", adata.var_names.is_unique)


# In[ ]:





# In[ ]:





# ## laber transfer Mono

# In[104]:


import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
import pyarrow.feather as feather


# In[105]:


# 使用to_string()打印完整内容
cell_state_counts = adata.obs['cell_type'].value_counts()
print(cell_state_counts.to_string())


# In[106]:


# 仅选择特定的细胞类型（
target_cell_types = ['Monocyte CD14+', 'Monocyte CD16+']
adata_subset = adata[adata.obs['cell_type'].isin(target_cell_types)].copy()


# In[107]:


# 读取表达矩阵
reference_expression = feather.read_feather('data/4_Mono_reference_expression.feather')
# 将 'gene_name' 列设置为行索引
reference_expression.set_index('gene_name', inplace=True)
# 读取元数据
reference_metadata= feather.read_feather('data/4_Mono_reference_metadata.feather')
# 将 'gene_name' 列设置为行索引
reference_metadata.set_index('cell_ID', inplace=True)
# 将表达矩阵转换为 AnnData 对象
reference = sc.AnnData(reference_expression.T)  # 转置矩阵以确保细胞为行，基因为列
reference


# In[108]:


# 将元数据添加到 AnnData 对象
reference.obs = reference_metadata

# 确保数据的维度对齐：细胞数（行数）应匹配元数据的行数
print(reference.obs.shape)  # 应该是 (27884, 84)
print(reference.shape)  # 应该是 (27884, 36601)
reference


# In[109]:


common_genes = reference.var_names.intersection(adata_subset.var_names)
# 截取 reference 的共同基因子集
reference = reference[:, common_genes].copy()
adata_subset = adata_subset[:, common_genes].copy()

assert all(reference.var_names == adata_subset.var_names), "Gene names are not aligned!"


# In[110]:


import numpy as np
print(np.any(reference.X < 0))  # 检查是否有负值
print(np.any(np.isnan(reference.X)))  # 检查是否有NaN


# In[ ]:


# 对两者都执行标准化、对数转换和 PCA
sc.pp.normalize_total(reference, target_sum=1e4)
sc.pp.log1p(reference)
sc.pp.scale(reference)
sc.tl.pca(reference)

sc.pp.normalize_total(adata_subset, target_sum=1e4)
sc.pp.log1p(adata_subset)
sc.pp.scale(adata_subset)
sc.tl.pca(adata_subset)


# In[113]:


# 运行 PCA 并构建邻接图
sc.pp.neighbors(reference, n_neighbors=15, use_rep='X_pca')
#sc.tl.umap(reference)  # 计算 UMAP 嵌入
sc.tl.umap(reference)
# 执行标签转移
sc.tl.ingest(adata_subset, reference, obs='ct_level_4')

# 将预测标签存储到查询数据
adata_subset.obs['predicted.id'] = adata_subset.obs['ct_level_4']

# 导出结果
adata_subset.obs.to_csv('data/4_Mono_Lable_transfer_predicted_neighbors.csv')


# ### PCA转化 trans_laber

# In[95]:


# 对参考数据和查询数据进行标准化
sc.pp.normalize_total(reference, target_sum=1e4)
sc.pp.log1p(reference)
#sc.pp.highly_variable_genes(reference)
sc.pp.highly_variable_genes(reference,  n_top_genes=2000)
sc.pp.scale(reference)

sc.pp.normalize_total(adata_subset, target_sum=1e4)
sc.pp.log1p(adata_subset)
#sc.pp.highly_variable_genes(adata_subset)
sc.pp.highly_variable_genes(adata_subset,  n_top_genes=2000)
sc.pp.scale(adata_subset)


# In[96]:


# 执行PCA降维
sc.tl.pca(reference, n_comps=50)
sc.tl.pca(adata_subset, n_comps=50)


# 提取PCA的结果
reference_pca = reference.obsm['X_pca']
query_pca = adata_subset.obsm['X_pca']

# 使用 NearestNeighbors 找到最近的锚点
nn = NearestNeighbors(n_neighbors=5)
nn.fit(reference_pca)  # 基于参考数据进行拟合
distances, indices = nn.kneighbors(query_pca)  # 查找查询数据中的最近邻

# 获取锚点 (anchors)
anchors = {
    'indices': indices,
    'distances': distances
}

# 获取参考数据的细胞类型标签
reference_labels = reference.obs['ct_level_4'].values

# 根据最近邻锚点进行标签转移
predicted_labels = [reference_labels[anchor_indices] for anchor_indices in indices]

import numpy as np
from sklearn.preprocessing import LabelEncoder

# 假设 predicted_labels 是一个列表，包含每个查询样本的预测标签（可能是字符串标签）
# 首先将参考数据的标签（如 reference_labels）传递给 LabelEncoder 并进行拟合
label_encoder = LabelEncoder()

# 假设 reference_labels 包含参考数据的所有可能的标签
label_encoder.fit(reference_labels)  # 在参考数据上拟合标签

# 将 predicted_labels 转换为整数标签
predicted_labels_encoded = [label_encoder.transform(labels) for labels in predicted_labels]

# 使用投票法选择每个查询样本的预测标签
predicted_ids = [np.bincount(labels).argmax() for labels in predicted_labels_encoded]

# 将预测标签添加到查询数据的 metadata
adata_subset.obs['predicted.id'] = label_encoder.inverse_transform(predicted_ids)


# In[97]:


adata_subset.obs.to_csv('data/4_Mono_Lable_transfer_predicted.csv')          


# ### figures

# In[98]:


# 打印预测结果的前几行
print(adata_subset.obs[['predicted.id']].head())

# 查看所有预测标签的分布
print(adata_subset.obs['predicted.id'].value_counts())

# 可视化查看预测标签的分布（如果你使用 seaborn 或 matplotlib）
import seaborn as sns
import matplotlib.pyplot as plt

sns.countplot(data=adata_subset.obs, x='predicted.id')
plt.xticks(rotation=45)
plt.title("Distribution of Predicted Labels")
plt.show()


# In[99]:


import seaborn as sns
import matplotlib.pyplot as plt

# 绘制predicted.id与time_point的分布图
plt.figure(figsize=(10, 6))
sns.countplot(data=adata_subset.obs, x='predicted.id', hue='time_point')

# 设置图表标题和标签
plt.title('Distribution of Predicted Labels by Time Point')
plt.xlabel('Predicted Label')
plt.ylabel('Count')

# 显示图形
plt.xticks(rotation=45)  # 旋转x轴标签以便查看
plt.legend(title='Time Point')
plt.show()


# In[ ]:





# In[ ]:





# ## laber transfer CD8T

# In[19]:


import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
import pyarrow.feather as feather


# In[15]:


# 打印出 'cell_state' 每个种类的数量
cell_state_counts = adata.obs['cell_state'].value_counts()
print(cell_state_counts)


# In[17]:


# 仅选择特定的细胞类型（T CD8 Naive, T CD8 Memory, T CD8 CTL）
target_cell_types = ['T CD8 Naive', 'T CD8 Memory', 'T CD8 CTL']
adata_subset = adata[adata.obs['cell_type'].isin(target_cell_types)].copy()


# In[30]:


# 读取表达矩阵
reference_expression = feather.read_feather('data/4_CD8T_reference_expression.feather')
# 将 'gene_name' 列设置为行索引
reference_expression.set_index('gene_name', inplace=True)
# 读取元数据
reference_metadata= feather.read_feather('data/4_CD8T_reference_metadata.feather')
# 将 'gene_name' 列设置为行索引
reference_metadata.set_index('cell_ID', inplace=True)


# In[33]:


# 检查数据维度
print(reference_metadata.shape)  # 检查元数据行数
print(reference_expression.shape)  # 检查表达矩阵行数


# In[34]:


import scanpy as sc

# 将表达矩阵转换为 AnnData 对象
reference = sc.AnnData(reference_expression.T)  # 转置矩阵以确保细胞为行，基因为列

# 将元数据添加到 AnnData 对象
reference.obs = reference_metadata

# 确保数据的维度对齐：细胞数（行数）应匹配元数据的行数
print(reference.obs.shape)  # 应该是 (27884, 84)
print(reference.shape)  # 应该是 (27884, 36601)


# In[ ]:


common_genes = reference.var_names.intersection(adata_subset.var_names)
# 截取 reference 的共同基因子集
reference = reference[:, common_genes].copy()
adata_subset = adata_subset[:, common_genes].copy()

assert all(reference.var_names == adata_subset.var_names), "Gene names are not aligned!"


# In[35]:


# 对参考数据和查询数据进行标准化
sc.pp.normalize_total(reference, target_sum=1e4)
sc.pp.log1p(reference)
#sc.pp.highly_variable_genes(reference)
sc.pp.highly_variable_genes(reference,  n_top_genes=2000)
sc.pp.scale(reference)

sc.pp.normalize_total(adata_subset, target_sum=1e4)
sc.pp.log1p(adata_subset)
#sc.pp.highly_variable_genes(adata_subset)
sc.pp.highly_variable_genes(adata_subset,  n_top_genes=2000)
sc.pp.scale(adata_subset)

# 执行PCA降维
sc.tl.pca(reference, n_comps=50)
sc.tl.pca(adata_subset, n_comps=50)


# In[36]:


# 提取PCA的结果
reference_pca = reference.obsm['X_pca']
query_pca = adata_subset.obsm['X_pca']

# 使用 NearestNeighbors 找到最近的锚点
nn = NearestNeighbors(n_neighbors=5)
nn.fit(reference_pca)  # 基于参考数据进行拟合
distances, indices = nn.kneighbors(query_pca)  # 查找查询数据中的最近邻

# 获取锚点 (anchors)
anchors = {
    'indices': indices,
    'distances': distances
}


# In[42]:


adata_subset


# In[53]:


adata_subset.var['highly_variable'].head()


# In[52]:


reference.var['highly_variable'].head()


# In[41]:


reference


# In[ ]:


# 获取参考数据的细胞类型标签
reference_labels = reference.obs['ct_level_4'].values

# 根据最近邻锚点进行标签转移
predicted_labels = [reference_labels[anchor_indices] for anchor_indices in indices]


# In[43]:


import numpy as np
from sklearn.preprocessing import LabelEncoder

# 假设 predicted_labels 是一个列表，包含每个查询样本的预测标签（可能是字符串标签）
# 首先将参考数据的标签（如 reference_labels）传递给 LabelEncoder 并进行拟合
label_encoder = LabelEncoder()

# 假设 reference_labels 包含参考数据的所有可能的标签
label_encoder.fit(reference_labels)  # 在参考数据上拟合标签

# 将 predicted_labels 转换为整数标签
predicted_labels_encoded = [label_encoder.transform(labels) for labels in predicted_labels]

# 使用投票法选择每个查询样本的预测标签
predicted_ids = [np.bincount(labels).argmax() for labels in predicted_labels_encoded]

# 将预测标签添加到查询数据的 metadata
adata_subset.obs['predicted.id'] = label_encoder.inverse_transform(predicted_ids)


# In[44]:


# 打印预测结果的前几行
print(adata_subset.obs[['predicted.id']].head())

# 查看所有预测标签的分布
print(adata_subset.obs['predicted.id'].value_counts())

# 可视化查看预测标签的分布（如果你使用 seaborn 或 matplotlib）
import seaborn as sns
import matplotlib.pyplot as plt

sns.countplot(data=adata_subset.obs, x='predicted.id')
plt.xticks(rotation=45)
plt.title("Distribution of Predicted Labels")
plt.show()


# In[45]:


adata_subset.obs.to_csv('data/4_CD8T_Lable_transfer_predicted.csv')


# In[46]:


import seaborn as sns
import matplotlib.pyplot as plt

# 绘制predicted.id与time_point的分布图
plt.figure(figsize=(10, 6))
sns.countplot(data=adata_subset.obs, x='predicted.id', hue='time_point')

# 设置图表标题和标签
plt.title('Distribution of Predicted Labels by Time Point')
plt.xlabel('Predicted Label')
plt.ylabel('Count')

# 显示图形
plt.xticks(rotation=45)  # 旋转x轴标签以便查看
plt.legend(title='Time Point')
plt.show()


# In[ ]:





# ## cell2tcr

# In[6]:


import cell2tcr
import celltypist
import scanpy as sc
import pandas as pd


# In[4]:


unique_values = adata.obs['cell_state'].unique()
print(unique_values)


# In[5]:


import scanpy as sc

# 归一化数据 (每个细胞的总计数归一化到 10000)
sc.pp.normalize_total(adata, target_sum=1e4)

# log1p 转换
sc.pp.log1p(adata)


# In[7]:


from celltypist import models
predictions = celltypist.annotate(adata, model='COVID19_HumanChallenge_Blood.pkl')
#COVID19_Immune_Landscape.pkl
#Healthy_COVID19_PBMC.pkl
adata = predictions.to_adata()
#adata = adata[adata.obs.predicted_labels.str.contains('T ')]
adata


# In[8]:


#标记激活的细胞
adata.obs['activated'] = adata.obs.predicted_labels.str.contains('Activated') & ~ adata.obs.predicted_labels.str.contains('MAI')

# 可视化结果
adata[adata.obs.activated].obs.predicted_labels.value_counts().plot.barh(title='Activated cells', xlabel='Count')


# In[9]:


import matplotlib.pyplot as plt

# 生成条形图
plt.figure(figsize=(8, 6))  # 设置图的大小，单位是英寸
adata[adata.obs.activated].obs.predicted_labels.value_counts().plot.barh()
plt.title('Activated All_data cells')  # 设置标题
plt.xlabel('Count')  # 设置x轴标签

# 保存为PDF文件
plt.savefig('figures/Tcell_data_activated_cells.pdf', bbox_inches='tight')


# In[15]:


import pandas as pd

# 将数据转换为 DataFrame
df = pd.DataFrame(adata.obs)

# 保存为 CSV 文件
df.to_csv('data/Tcell_data_activated_cells_Challenge.csv', index=False)


# In[ ]:





# ## Umap 图 绘制

# In[11]:


import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt


# In[12]:


print(adata.obs['cell_compartment'].unique())


# In[13]:


# 确保 ct_merge 是字符串类型，避免与 Categorical 类型冲突
adata.obs['ct_merge'] = adata.obs['cell_compartment'].astype(str)

# 修改细胞类型
adata.obs['ct_merge'] = adata.obs['ct_merge'].replace({
    'T CD4': 'CD4T',
    'T Reg': 'CD4T',
    'T CD8': 'CD8T',
    'T G/D': 'gdT',
    'T MAI': 'MAIT',
    'T DN': 'Tdn',
    'NK': 'NK',
    'NK CD56+': 'NK',
    'B': 'Bcell',
    'B ABS': 'Bcell',
    'HPC': 'HPSC',
})

# 使用 loc 方法进一步修改特定条件的值
adata.obs.loc[adata.obs['cell_type'] == 'pDC', 'ct_merge'] = 'pDC'
adata.obs.loc[adata.obs['cell_type'].isin(['cDC1', 'cDC2', 'cDC3', 'AS-DC']), 'ct_merge'] = 'cDC'

# 将 ct_merge 转回 Categorical 类型（如果需要）
adata.obs['ct_merge'] = pd.Categorical(adata.obs['ct_merge'])

# 验证修改后的独特细胞类型
print(adata.obs['ct_merge'].unique())


# In[14]:


# 假设 adata 是 AnnData 对象
import pandas as pd

# 筛选 'CD4T' 和 'CD8T' 的总数量
total_count = adata.obs['ct_merge'].isin(['CD4T', 'CD8T']).sum()

# 打印统计结果
print(f"Total count of CD4T and CD8T: {total_count}")


# In[8]:


# 设置时间点和COVID状态的因子顺序（可选）
adata.obs['time_point'] = pd.Categorical(
    adata.obs['time_point'], categories=['D-1', 'D3', 'D7', 'D10', 'D14', 'D28'], ordered=True
)
adata.obs['covid_status'] = pd.Categorical(
    adata.obs['covid_status'], categories=['Sustained infection', 'Transient infection', 'Abortive infection'], ordered=True
)


# In[31]:


import matplotlib as mpl
import matplotlib.pyplot as plt

# 设置字体为 sans-serif 类型
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']  # R 通常用 Helvetica 或类似字体
mpl.rcParams['pdf.fonttype'] = 42  # 确保 PDF 中嵌入矢量字体

def plot_umap_with_labels(adata, color_col, palette, output_file, label_size=12, width=28, height=20):
    """
    绘制UMAP图，两种布局：基础UMAP图和带标签UMAP图。
    
    Parameters:
        adata: AnnData 对象
        color_col: str, 用于着色的列
        palette: dict, 颜色映射
        output_file: str, 保存图像的路径（PDF 文件）
        label_size: int, 标签字体大小
        width: int, 图像宽度
        height: int, 图像高度
    """
    if 'X_umap_harmony_rna_wvdj_30pcs_6000hvgs' not in adata.obsm.keys():
        raise ValueError("UMAP embedding not found in .obsm. Please ensure it is computed.")
    
    umap_df = pd.DataFrame(
        adata.obsm['X_umap_harmony_rna_wvdj_30pcs_6000hvgs'],
        columns=['umap_1', 'umap_2'],
        index=adata.obs.index
    )
    umap_df[color_col] = adata.obs[color_col]

    cell_type_med = umap_df.groupby(color_col)[['umap_1', 'umap_2']].median().reset_index()

    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(output_file) as pdf:
        # 绘制基础UMAP图
        plt.figure(figsize=(width, height))
        unique_categories = umap_df[color_col].unique()
        for cat in unique_categories:
            subset = umap_df[umap_df[color_col] == cat]
            plt.scatter(
                subset['umap_1'], subset['umap_2'], 
                label=cat, s=5, alpha=0.8, color=palette.get(cat, 'grey')
            )
        
        plt.title(f'UMAP colored by {color_col}')
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.legend(markerscale=2, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        pdf.savefig()
        plt.close()

        # 绘制带标签的UMAP图
        plt.figure(figsize=(width, height))
        texts = []
        for cat in unique_categories:
            subset = umap_df[umap_df[color_col] == cat]
            plt.scatter(
                subset['umap_1'], subset['umap_2'], 
                label=cat, s=5, alpha=0.8, color=palette.get(cat, 'grey')
            )
        
        for _, row in cell_type_med.iterrows():
            texts.append(
                plt.text(
                    row['umap_1'], row['umap_2'], 
                    row[color_col], fontsize=label_size, fontweight='bold',
                    color=palette.get(row[color_col], 'grey'),
                    bbox=dict(boxstyle="round,pad=0.3", edgecolor=palette.get(row[color_col], 'grey'), facecolor='white')
                )
            )
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', lw=0.5))
        plt.title(f'UMAP with {color_col} labels')
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.tight_layout()
        pdf.savefig()
        plt.close()


# In[33]:


# 调用函数
colpalette = {
    'Bcell': '#AB3282',
    'CD4T': '#6778AE',
    'CD8T': '#53A85F',
    'cDC': '#9FA3A8',
    'gdT': '#23452F',
    'HPSC': '#585658',
    'ILC': '#F1BC62',
    'MAIT': '#3A6963',
    'Monocyte': '#E95C39',
    'NK': '#E1A111',
    'pDC': '#938175',
    'Tdn': '#495608'
}

plot_umap_with_labels(
    adata, 
    color_col='ct_merge', 
    palette=colpalette, 
    output_file='0_umap_challenge_ct_merge.pdf',
    label_size=72
)


# ##  ATG7基因

# In[3]:


ATG7= ["ATG7"]
sc.tl.score_genes(adata, ATG7, 
                  ctrl_size=50, gene_pool=None,n_bins=25,
                  score_name='ATG71', random_state=0, copy=False, use_raw=None)


# In[ ]:





# In[31]:


# import pandas as pd

# # 读取 CSV 文件，将行名作为普通列
# df_genes = pd.read_csv(
#     "/share/home/qlab/projects/project_cvv/yyp_results_3/5_Cell_type_level_3/2_Mono/7_diff_gene_in_Mono/data/macroautophagy_genes.csv",
#     index_col=0
# )

# # 重置行名为普通列，并命名为 'X'
# df_genes = df_genes.reset_index().rename(columns={"index": "X"})
# macroautophagy_genes = df_genes['X']

# 查看基因数量
# print("Number of genes:", len(macroautophagy_genes))


# In[11]:


macroautophagy_genes =["ATG3", "ATG10", "ATG12", "ATG5", "ATG4",
                         "MAP1LC3B", "ATG16L1", "BECN1", "ULK1"]
# ["ATG7","ATG3", "ATG10", "ATG12", "ATG5", "ATG4", "MAP1LC3B",
#                           "ATG16L1", "BECN1", "ULK1", "ATG9A", "ATG9B", "WIPI1",
#                           "WIPI2", "GABARAP", "FIP200", "PIK3C3", "SQSTM1",  
#                           "NBR1", "ATG14", "TECPR1", "BNIP3"]

sc.tl.score_genes(adata, macroautophagy_genes, 
                  ctrl_size=50, gene_pool=None,n_bins=25,
                  score_name='mac_ATG71', random_state=0, copy=False, use_raw=None)


# ## pDC 凋亡基因

# In[5]:


genes_Apoptosis = ["AGO3", "CBX4", "UBC", "DAPK1", "TNFRSF1B",
                      "MDM4", "PIK3R1", "FOS", "SIVA1", "APPL1", 
                      "RAF1", "JAK2", "HIST1H2AC", "HIST1H2BC", 
                      "HIPK1", "HIST1H2BK", "CDKN2D", "JUN", "CBX6", 
                      "PHC3", "MOAP1", "BAK1", "GSK3B", "RPS27A"]
sc.tl.score_genes(adata, genes_Apoptosis, 
                  ctrl_size=50, gene_pool=None,n_bins=25,
                  score_name='Apoptosis1', random_state=0, copy=False, use_raw=None)


# In[ ]:





# ## VIA gene Module

# In[6]:


genes_VIA = ["ACTG1", "PCLAF", "ACTB", "GAPDH", "CD38", "GZMH", "PFN1", "CORO1A", "PYCARD", "PSMB8", "TMSB10", "ARPC3", 
  "CLIC1", "CD52", "TMSB4X", "APOBEC3H", "ADA", "SEM1", "MYL6", "GZMA", "DECR1", "ARPC2", "NOP10", "APOBEC3C", 
  "LGALS1", "GZMB", "CFL1", "CLTA", "TWF2", "ISG15", "MXD4", "SUB1", "PSMA4", "HLA-DRA", "PPP1R18", "GSTO1", 
  "APOBEC3G", "ARPC5", "UCP2", "PSMA5", "WDR1", "ARHGDIB", "YWHAE", "CCL5", "SLC25A5", "PSMB4", "ITGB1", 
  "DYNLL1", "ARPC4", "ATP5MC2", "FGFBP2", "GTF3C6", "RGS10", "ATP6V1D", "S1PR4", "EZH2", "ATP5F1A", "GMFG", 
  "OSTF1", "PRF1", "ABRACL", "GLRX", "ACTR3", "SIT1", "TXNDC17", "H3F3A", "CCDC28B", "ATP5MG", "DNAJC15", 
  "IDH2", "SEC11A", "FABP5", "CAPZB", "DBI", "IFI16", "DCTN2", "COTL1", "CNN2", "CHMP2A", "ATP5MC3", "ARPC1B", 
  "NDUFS2", "MT1E", "MRPL28", "OAS2", "SERPINB1", "CARHSP1", "TIMD4", "UBE2L6", "MYL6B", "DCTN3", "IRF4", 
  "ANXA5", "ATP5F1B", "GIMAP4", "PSMB2", "NDUFB3", "NCOA4", "YARS", "CARD16", "COX5A", "FKBP1A", "TCEAL8", 
  "CLDND1", "NDUFB5", "ATP5F1C", "MRPL42", "RPA3", "MT1F", "NUDT5", "POLR2G", "LDHB", "GBP1", "ATP5MF", 
  "COX6B1", "COX6C", "SUCLG1", "TAP1", "CLIC3", "TPM4", "GBP2", "NDUFC2", "CHCHD1", "PPP1CA", "TKT", "SH3BGRL3", 
  "EIF4E2", "LAIR2", "MRPL51", "TRAPPC1", "CSK", "BCL2L11", "LIMD2", "JPT1", "GTF2H5", "NDUFA12", "TROAP", 
  "RTRAF", "ZYX", "RAC2", "CHCHD5", "MPG", "RPS6KA1", "SRSF9", "CEBPD", "GIMAP1", "COPZ1", "TALDO1", "CALM3", 
  "EIF2S2", "HLA-DMA", "POMP", "PTRHD1", "COMMD4", "GGCT", "HAVCR2", "LAP3", "FERMT3", "COPS9", "QARS", 
  "AP2S1", "ZNHIT1", "ICA1", "HADHB", "BLOC1S1", "RABL3", "IFI27L2", "RACK1", "CD27", "FKBP3", "MRPL10", 
  "TXNDC9", "FIBP", "PTTG1", "PSMD8", "CXCR3", "PSME1", "LCP1", "ACAA2", "CHI3L2", "NCKAP1L", "PPP4C", "IGBP1", 
  "PSMB3", "LAMTOR2", "ARHGAP30", "TMEM256", "COPB2", "AC010618.1", "AP1S1", "PFDN4", "NT5C3A", "HMGA1", 
  "MT-CO1", "PSMA2", "NSMCE1", "PARK7"]
sc.tl.score_genes(adata, genes_VIA, 
                  ctrl_size=50, gene_pool=None,n_bins=25,
                  score_name='VIA1', random_state=0, copy=False, use_raw=None)




# ## Tcell

# In[1]:


import anndata as ad

# 读取数据
adata = ad.read_h5ad("data/challenge_pbmc_cellxgene_230223.h5ad")

# 筛选出cell_compartment列中以"T"开头的细胞数据
Tcell_data = adata[adata.obs['cell_compartment'].str.startswith('T')].copy()

# 查看结果
print(Tcell_data)


# In[4]:


import pandas as pd
import scanpy as sc
antigen_module_genes_df = pd.read_csv("/share/home/qlab/projects/project_cvv/yyp_results_3/0_paper_data/2_VIA_paper_data/data/antigen_module_genes.csv")
antigen_module_genes = antigen_module_genes_df["x"].tolist()  # 假设 genes 列包含基因名称
sc.tl.score_genes(Tcell_data, antigen_module_genes, 
                  ctrl_size=50, gene_pool=None,n_bins=25,
                  score_name='VIAScore', random_state=0, copy=False, use_raw=None)


# In[6]:


# 保存为 CSV 文件
Tcell_data.obs.to_csv('data/Tcell_challenge_adata_obs.csv', index=True)


# In[97]:


# # 确定 time_point 的顺序
# time_point_order = ["D-1", "D3", "D7", "D10", "D14", "D28"]

# # 将 time_point 列转换为类别类型，并设置顺序
# adata.obs['time_point'] = pd.Categorical(adata.obs['time_point'], categories=time_point_order, ordered=True)

adata.obs['time_point'].value_counts()


# ## group line

# In[34]:


covid_status_values


# In[36]:


subset_df['covid_status'].unique()


# In[37]:


# 提取元数据
meta_data_df = adata.obs.copy()

# 确定 time_point 的顺序
time_point_order = ["D-1", "D3", "D7", "D10", "D14", "D28"]

# 将 time_point 列转换为类别类型，并设置顺序
meta_data_df['time_point'] = pd.Categorical(meta_data_df['time_point'], categories=time_point_order, ordered=True)

# 创建颜色调色板
colpalette_27 = sns.color_palette("Set1", 27)

# 获取所有 covid_status 的值
covid_status_values = meta_data_df['covid_status'].unique()

for status in covid_status_values:
    # 过滤数据
    subset_df = meta_data_df[meta_data_df['covid_status'] == status]
    # 计算 Celltype_prop
    subset_df['Count_time_point'] = subset_df.groupby('time_point', observed=False)['time_point'].transform('count')
    subset_df['Count_time_point_celltype'] = subset_df.groupby(['time_point', 'cell_type'], observed=False)['time_point'].transform('count')
    subset_df['Celltype_prop'] = subset_df['Count_time_point_celltype'] / subset_df['Count_time_point']
    # 绘制分面折线图
    g = sns.FacetGrid(subset_df, col="cell_type", col_wrap=6, height=3, aspect=1.5, sharey=False)
    g.map_dataframe(sns.lineplot, x='time_point', y='Celltype_prop', marker="o", hue='cell_type', palette=colpalette_27)

    # 添加图例和调整布局
    g.add_legend()
    g.set_axis_labels("Time Point", "Proportion")
    g.fig.suptitle(f"Line plot for time_point vs cell_type proportion\n(COVID Status: {status})", y=1.02)

    # 保存图表
    g.savefig(f"1_Challenge_Line_time_point_cell_type_{status.replace(' ', '_')}.pdf", bbox_inches='tight')
    plt.close()


# In[22]:


# 提取元数据
meta_data_df = adata.obs.copy()

# 创建颜色调色板
colpalette_27 = sns.color_palette("Set1",27)
colpalette_64 = sns.color_palette("Set1", 64)
# 确定 time_point 的顺序
time_point_order = ["D-1", "D3", "D7", "D10", "D14", "D28"]

# 将 time_point 列转换为类别类型，并设置顺序
meta_data_df['time_point'] = pd.Categorical(meta_data_df['time_point'], categories=time_point_order, ordered=True)


# In[ ]:


# 计算 Celltype_prop
meta_data_df['Count_time_point'] = meta_data_df.groupby('time_point', observed=False)['time_point'].transform('count')
meta_data_df['Count_time_point_celltype'] = meta_data_df.groupby(['time_point', 'cell_type'], observed=False)['time_point'].transform('count')
meta_data_df['Celltype_prop'] = meta_data_df['Count_time_point_celltype'] / meta_data_df['Count_time_point']

# 绘制分面折线图
g = sns.FacetGrid(meta_data_df, col="cell_type", col_wrap=6, height=3, aspect=1.5, sharey=False)
g.map_dataframe(sns.lineplot, x='time_point', y='Celltype_prop', marker="o", hue='cell_type', palette=colpalette_64)

# 添加图例和调整布局
g.add_legend()
g.set_axis_labels("Time Point", "Proportion")
g.fig.suptitle("Line plot for time_point vs cell_type proportion", y=1.02)

# 保存图表
g.savefig("1_Challenge_Line_time_point_cell_cell_type.pdf", bbox_inches='tight')
plt.close()


# In[26]:


# 计算 Celltype_prop
meta_data_df['Count_time_point'] = meta_data_df.groupby('time_point', observed=False)['time_point'].transform('count')
meta_data_df['Count_time_point_celltype'] = meta_data_df.groupby(['time_point', 'cell_state_woIFN'], observed=False)['time_point'].transform('count')
meta_data_df['Celltype_prop'] = meta_data_df['Count_time_point_celltype'] / meta_data_df['Count_time_point']
# 绘制分面折线图
g = sns.FacetGrid(meta_data_df, col="cell_type", col_wrap=6, height=3, aspect=1.5, sharey=False)
g.map_dataframe(sns.lineplot, x='time_point', y='Celltype_prop', marker="o", hue='cell_state_woIFN', palette=colpalette_64)
# 添加图例和调整布局
g.add_legend()
g.set_axis_labels("Time Point", "Proportion")
g.fig.suptitle("Line plot for time_point vs cell_type proportion", y=1.02)
# 保存图表
g.savefig("1_Challenge_Line_time_point_cell_state_woIFN.pdf", bbox_inches='tight')
plt.close()



# 计算 Celltype_prop
meta_data_df['Count_time_point'] = meta_data_df.groupby('time_point', observed=False)['time_point'].transform('count')
meta_data_df['Count_time_point_celltype'] = meta_data_df.groupby(['time_point', 'cell_state_woIFN'], observed=False)['time_point'].transform('count')
meta_data_df['Celltype_prop'] = meta_data_df['Count_time_point_celltype'] / meta_data_df['Count_time_point']
# 绘制分面折线图
g = sns.FacetGrid(meta_data_df, col="cell_state_woIFN", col_wrap=8, height=3, aspect=1.5, sharey=False)
g.map_dataframe(sns.lineplot, x='time_point', y='Celltype_prop', marker="o", hue='cell_state_woIFN', palette=colpalette_64)
# 添加图例和调整布局
g.add_legend()
g.set_axis_labels("Time Point", "Proportion")
g.fig.suptitle("Line plot for time_point vs cell_type proportion", y=1.02)
# 保存图表
g.savefig("1_Challenge_Line_time_point_cell_state_woIFN_2.pdf", bbox_inches='tight')
plt.close()


# In[29]:


colpalette_83 = sns.color_palette("Set1",83)

# 计算 Celltype_prop
meta_data_df['Count_time_point'] = meta_data_df.groupby('time_point', observed=False)['time_point'].transform('count')
meta_data_df['Count_time_point_celltype'] = meta_data_df.groupby(['time_point', 'cell_state'], observed=False)['time_point'].transform('count')
meta_data_df['Celltype_prop'] = meta_data_df['Count_time_point_celltype'] / meta_data_df['Count_time_point']
# 绘制分面折线图
g = sns.FacetGrid(meta_data_df, col="cell_type", col_wrap=6, height=3, aspect=1.5, sharey=False)
g.map_dataframe(sns.lineplot, x='time_point', y='Celltype_prop', marker="o", hue='cell_state', palette=colpalette_83)
# 添加图例和调整布局
g.add_legend()
g.set_axis_labels("Time Point", "Proportion")
g.fig.suptitle("Line plot for time_point vs cell_type proportion", y=1.02)
# 保存图表
g.savefig("1_Challenge_Line_time_point_cell_state_cell_state.pdf", bbox_inches='tight')
plt.close()



# 计算 Celltype_prop
meta_data_df['Count_time_point'] = meta_data_df.groupby('time_point', observed=False)['time_point'].transform('count')
meta_data_df['Count_time_point_celltype'] = meta_data_df.groupby(['time_point', 'cell_state'], observed=False)['time_point'].transform('count')
meta_data_df['Celltype_prop'] = meta_data_df['Count_time_point_celltype'] / meta_data_df['Count_time_point']
# 绘制分面折线图
g = sns.FacetGrid(meta_data_df, col="cell_state", col_wrap=8, height=3, aspect=1.5, sharey=False)
g.map_dataframe(sns.lineplot, x='time_point', y='Celltype_prop', marker="o", hue='cell_state', palette=colpalette_83)
# 添加图例和调整布局
g.add_legend()
g.set_axis_labels("Time Point", "Proportion")
g.fig.suptitle("Line plot for time_point vs cell_type proportion", y=1.02)
# 保存图表
g.savefig("1_Challenge_Line_time_point_cell_state_2.pdf", bbox_inches='tight')
plt.close()


# ## IFN HLA VIA module

# In[7]:


genes_ifn = ["IFI35", "IFI44L", "IFI6", "IFIT3", 
             "IRF7", "ISG15", "MX1", "MX2", 
             "OAS1", "OAS2"]
sc.tl.score_genes(adata, genes_ifn, 
                  ctrl_size=50, gene_pool=None,n_bins=25,
                  score_name='ifn_stim', random_state=0, copy=False, use_raw=None)


# In[8]:


HLA = ["HLA-DQA2"]
sc.tl.score_genes(adata, HLA, 
                  ctrl_size=50, gene_pool=None,n_bins=25,
                  score_name='HLA_DQA2', random_state=0, copy=False, use_raw=None)


# In[9]:


adata


# In[12]:


# 保存为 CSV 文件
adata.obs.to_csv('data/challenge_adata_obs.csv', index=True)


# # figures

# In[104]:


# 创建画布和子图
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# 获取covid_status的所有唯一值
statuses = adata.obs['covid_status'].unique()

# 遍历每个covid_status并在对应的子图上绘制
for i, status in enumerate(statuses):
    ax = axes.flatten()[i]
    sc.pl.umap(adata[adata.obs['covid_status'] == status], vmax= 20,vmin= -20 ,color='VIA_score', cmap='coolwarm', ax=ax, show=False)
    ax.set_title(f'COVID Status: {status}')

# 保存图像
plt.savefig('2_Challenge_Umap_VIA_score_covid_status.pdf')

# 显示图像
plt.show()


# In[103]:


# 创建画布和子图
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# 获取covid_status的所有唯一值

statuses = ["D-1", "D3", "D7", "D10", "D14", "D28"]

# 遍历每个covid_status并在对应的子图上绘制
for i, status in enumerate(statuses):
    ax = axes.flatten()[i]
    sc.pl.umap(adata[adata.obs['time_point'] == status],vmax= 20,vmin= -20 , color='VIA_score', cmap='coolwarm', ax=ax, show=False)
    ax.set_title(f'COVID Status: {status}')

# 保存图像
plt.savefig('2_Challenge_Umap_VIA_score_time_point.pdf')

# 显示图像
plt.show()


# In[105]:


# 创建画布和子图
fig, axes = plt.subplots(12, 7, figsize=(50, 35))

# 获取covid_status的所有唯一值
statuses = adata.obs['cell_state'].unique()

# 遍历每个covid_status并在对应的子图上绘制
for i, status in enumerate(statuses):
    ax = axes.flatten()[i]
    sc.pl.umap(adata[adata.obs['cell_state'] == status], vmax= 20,vmin= -20 ,color='VIA_score', cmap='coolwarm', ax=ax, show=False)
    ax.set_title(f'COVID Status: {status}')

# 保存图像
plt.savefig('2_Challenge_Umap_VIA_score_cell_state.pdf')

# 显示图像
plt.show()


# In[102]:


# 创建画布和子图
fig, axes = plt.subplots(9, 3, figsize=(15, 45))

# 获取covid_status的所有唯一值
statuses = adata.obs['cell_type'].unique()

# 遍历每个covid_status并在对应的子图上绘制
for i, status in enumerate(statuses):
    ax = axes.flatten()[i]
    sc.pl.umap(adata[adata.obs['cell_type'] == status], vmax= 20,vmin= -20 ,color='VIA_score', cmap='coolwarm', ax=ax, show=False)
    ax.set_title(f'COVID Status: {status}')

# 保存图像
plt.savefig('2_Challenge_Umap_VIA_score_cell_type.pdf')

# 显示图像
plt.show()


# In[101]:


# 设置 UMAP 坐标
adata.obsm['X_umap'] = adata.obsm['X_umap_harmony_rna_wvdj_30pcs_6000hvgs']

# 绘制 UMAP 图，使用 coolwarm 颜色映射
sc.pl.umap(adata, color='VIA_score', cmap='coolwarm',vmax= 20,vmin= -20 ,save='2_Challenge_Umap_VIA_score.pdf') #groups="covid_status",


# In[ ]:


# 提取元数据
meta_data_df = adata.obs.copy()

# 获取所有 covid_status 的值
covid_status_values = meta_data_df['covid_status'].unique()

for status in covid_status_values:
    # 过滤数据
    subset_df = meta_data_df[meta_data_df['covid_status'] == status]
    # 计算 Celltype_prop
    subset_df['Count_time_point'] = subset_df.groupby('time_point', observed=False)['time_point'].transform('count')
    subset_df['Count_time_point_celltype'] = subset_df.groupby(['time_point', 'cell_type'], observed=False)['time_point'].transform('count')
    subset_df['Celltype_prop'] = subset_df['Count_time_point_celltype'] / subset_df['Count_time_point']
    # 绘制分面折线图
    g = sns.FacetGrid(subset_df, col="cell_type", col_wrap=6, height=3, aspect=1.5, sharey=False)
    g.map_dataframe(sns.lineplot, x='time_point', y='Celltype_prop', marker="o", hue='cell_type', palette=colpalette_27)

    # 添加图例和调整布局
    g.add_legend()
    g.set_axis_labels("Time Point", "Proportion")
    g.fig.suptitle(f"Line plot for time_point vs cell_type proportion\n(COVID Status: {status})", y=1.02)

    # 保存图表
    g.savefig(f"1_Challenge_Line_time_point_cell_type_{status.replace(' ', '_')}.pdf", bbox_inches='tight')
    plt.close()


# ### ATG7 NELL2

# In[ ]:


c("ATG7",'NELL2')


# In[13]:


import scanpy as sc
import matplotlib.pyplot as plt

# 假设你的adata对象已经加载
# adata = sc.read_h5ad('你的adata文件.h5ad')

# 设置 UMAP 坐标
adata.obsm['X_umap'] = adata.obsm['X_umap_harmony_rna_wvdj_30pcs_6000hvgs']

# 绘制UMAP图并突出显示ATG7基因表达
sc.pl.umap(adata, color='ATG7',vmax= 3,vmin= -3 ,
           cmap='coolwarm',title='UMAP of ATG7 Expression',save='3_Challenge_ATG7.pdf')

# 绘制UMAP图并突出显示NELL2基因表达
sc.pl.umap(adata, color='NELL2',vmax= 3,vmin= -3 ,
           cmap='coolwarm',title='UMAP of NELL2 Expression',save='2_Challenge_NELL2.pdf')


#plt.show()


# In[14]:


# 创建画布和子图
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# 获取covid_status的所有唯一值
statuses = adata.obs['covid_status'].unique()

# 遍历每个covid_status并在对应的子图上绘制
for i, status in enumerate(statuses):
    ax = axes.flatten()[i]
    sc.pl.umap(adata[adata.obs['covid_status'] == status], vmax= 3,vmin= -3 ,color='ATG7', cmap='coolwarm', ax=ax, show=False)
    ax.set_title(f'COVID Status: {status}')

# 保存图像
plt.savefig('3_Challenge_ATG7.pdf')

# 显示图像
plt.show()


# In[15]:


# 创建画布和子图
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# 获取covid_status的所有唯一值
statuses = adata.obs['covid_status'].unique()

# 遍历每个covid_status并在对应的子图上绘制
for i, status in enumerate(statuses):
    ax = axes.flatten()[i]
    sc.pl.umap(adata[adata.obs['covid_status'] == status], vmax= 3,vmin= -3 ,color='NELL2', cmap='coolwarm', ax=ax, show=False)
    ax.set_title(f'COVID Status: {status}')

# 保存图像
plt.savefig('3_Challenge_NELL2.pdf')

# 显示图像
plt.show()


# In[ ]:




