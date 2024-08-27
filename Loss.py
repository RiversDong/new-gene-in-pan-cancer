'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import networkx as nx
from scipy.stats import pearsonr
import matplotlib as mpl

# svg可以编辑文字

mpl.rcParams["pdf.fonttype"]=42
#mpl.rcParams['font.family'] = 'sans-serif'
#mpl.rcParams['font.family'] = 'Arial'



def get_age(data, age_type):
    gene_age = data[["Gene", age_type]]
    gene_age = gene_age.dropna(axis=0, how="any")
    return gene_age

tcga = "/data/chuand/essential_gene_evo/data/TCGA_interaction_type_top_tcga_mutual_20"
gtex = "/data/chuand/essential_gene_evo/data/GTEx_interaction_type_top_gtex_mutual_20"
age = "/data/chuand/essential_gene_evo/data/13059_2022_2821_MOESM2_ESM_Ma_gb_2022.txt"
tcga = pd.read_csv(tcga, sep="\t")
gtex = pd.read_csv(gtex, sep="\t")
age_data = pd.read_csv(age, sep="\t", low_memory=False)
gene2age = get_age(age_data, "4 merged age groups")
gene2age = gene2age.set_index('Gene')['4 merged age groups'].to_dict()

def add_age(df, age_dict):
    add_age = []
    for index, row in df.iterrows():
        gene1, gene2 = row["gene1"], row["gene2"]
        age1, age2 = age_dict.get(gene1, "unknown"), age_dict.get(gene2, "unknown")
        add_age.append([gene1, gene2, age1, age2])
    new_coexp = pd.DataFrame(add_age, columns=["Gene1", "Gene2", "Age1", "Age2"])
    return new_coexp

new_gtex = add_age(gtex, gene2age)

gtex_list = [i + "_" + j for i, j in zip(gtex["gene1"], gtex["gene2"])]
tcga_list = [i + "_" + j for i, j in zip(tcga["gene1"], tcga["gene2"])]

out_res = "/data/chuand/essential_gene_evo/result/lossed.csv"
age_groups = ["UC", "EM", "MM", "PSG"]
infoes = pd.read_csv(out_res, sep=",")
head_map_data = []
gene2time = {"PSG": 73, "MM": 180, "EM": 758, "UC": 4250}
time_gaps = []
mean_ints = []
labels = []

within_age_group = {}

within_age_group["UC_VS_UC"] = 1746
within_age_group["EM_VS_EM"] = 289
within_age_group["MM_VS_MM"] = 53.5
within_age_group["PSG_VS_PSG"] = 36.5

for i in age_groups:
    i_data = []
    G = nx.Graph()
    for j in age_groups:
    
        if f"{i}_VS_{j}" in within_age_group:
            ij_data = infoes[((infoes["age1"] == i) & (infoes["age2"] == j)) | ((infoes["age1"] == j) & (infoes["age2"] == i))]
            edgeLen = len(ij_data)
            ij_data_gtex = new_gtex[((new_gtex["Age1"] == i) & (new_gtex["Age2"] == j)) | ((new_gtex["Age1"] == j) & (new_gtex["Age2"] == i))]
            all_edgeLen = len(ij_data_gtex)
            mean_int = round(edgeLen / all_edgeLen, 5) if all_edgeLen != 0 else 0
            i_data.append(mean_int)
            labels.append(f"{i}-{j}")
            time_gap = within_age_group[f"{i}_VS_{j}"]
            time_gaps.append(time_gap)
            mean_ints.append(mean_int)
            print(time_gap, mean_int, f"{i}_VS_{j}")  
        else:
            ij_data = infoes[((infoes["age1"] == i) & (infoes["age2"] == j)) | ((infoes["age1"] == j) & (infoes["age2"] == i))]
            edgeLen = len(ij_data)
            ij_data_gtex = new_gtex[((new_gtex["Age1"] == i) & (new_gtex["Age2"] == j)) | ((new_gtex["Age1"] == j) & (new_gtex["Age2"] == i))]
            all_edgeLen = len(ij_data_gtex)
            mean_int = round(edgeLen / all_edgeLen, 5) if all_edgeLen != 0 else 0
            i_data.append(mean_int)
            labels.append(f"{i}-{j}")
            time_gap = abs(gene2time[i] - gene2time[j])
            time_gaps.append(time_gap)
            mean_ints.append(mean_int)
            print(time_gap, mean_int, f"{i}_VS_{j}")
    head_map_data.append(i_data)

# 绘制 time_gaps 和 mean_ints 的散点图，并拟合直线，计算皮尔逊相关系数
plt.figure()

colors = plt.cm.viridis(np.linspace(0, 1, 16))

# 增加透明度
plt.scatter(time_gaps, mean_ints, color=colors, alpha=0.7)

for i, (time_gap, mean_int, label, color) in enumerate(zip(time_gaps, mean_ints, labels, colors)):
    plt.scatter(time_gap, mean_int, color=color, alpha=0.7, label=label)

plt.xlabel('Time gaps')
plt.ylabel('Loss percentage')
#plt.title('Time Gaps vs Mean Interaction Loss')

# 拟合直线
m, b = np.polyfit(time_gaps, mean_ints, 1)
plt.plot(time_gaps, m*np.array(time_gaps) + b, color='red')

# 计算皮尔逊相关系数
corr, p_value = pearsonr(time_gaps, mean_ints)
print(corr, p_value)
plt.text(max(time_gaps)*0.7, max(mean_ints)*0.9, f'R = {corr:.2f}')
plt.legend()
plt.tight_layout()
plt.savefig("/data/chuand/essential_gene_evo/result/time_gap_vs_mean_interaction_loss.pdf")
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from scipy.stats import pearsonr
import matplotlib as mpl

# SVG can edit text
mpl.rcParams["pdf.fonttype"] = 42

def get_age(data, age_type):
    gene_age = data[["Gene", age_type]]
    gene_age = gene_age.dropna(axis=0, how="any")
    return gene_age

# Load data
tcga = pd.read_csv("/data/chuand/essential_gene_evo/data/TCGA_interaction_type_top_tcga_mutual_20", sep="\t")
gtex = pd.read_csv("/data/chuand/essential_gene_evo/data/GTEx_interaction_type_top_gtex_mutual_20", sep="\t")
age_data = pd.read_csv("/data/chuand/essential_gene_evo/data/13059_2022_2821_MOESM2_ESM_Ma_gb_2022.txt", sep="\t", low_memory=False)

gene2age = get_age(age_data, "4 merged age groups")
gene2age = gene2age.set_index('Gene')['4 merged age groups'].to_dict()

def add_age(df, age_dict):
    add_age = []
    for index, row in df.iterrows():
        gene1, gene2 = row["gene1"], row["gene2"]
        age1, age2 = age_dict.get(gene1, "unknown"), age_dict.get(gene2, "unknown")
        add_age.append([gene1, gene2, age1, age2])
    new_coexp = pd.DataFrame(add_age, columns=["Gene1", "Gene2", "Age1", "Age2"])
    return new_coexp

new_gtex = add_age(gtex, gene2age)

out_res = "/data/chuand/essential_gene_evo/result/lossed.csv"
age_groups = ["UC", "EM", "MM", "PSG"]
infoes = pd.read_csv(out_res, sep=",")
gene2time = {"PSG": 73, "MM": 180, "EM": 758, "UC": 4250}
time_gaps = []
mean_ints = []
labels = []

within_age_group = {
    "UC_VS_UC": 1746,
    "EM_VS_EM": 289,
    "MM_VS_MM": 53.5,
    "PSG_VS_PSG": 36.5
}

# Combine symmetric pairs like UC-EM and EM-UC
processed_pairs = set()

for i in age_groups:
    for j in age_groups:
        pair = tuple(sorted([i, j]))  # Sort to handle symmetry, e.g., (UC, EM) == (EM, UC)
        if pair in processed_pairs:
            continue
        
        processed_pairs.add(pair)
        label = f"{pair[0]}-{pair[1]}({pair[1]}-{pair[0]})"
        
        if f"{i}_VS_{j}" in within_age_group or f"{j}_VS_{i}" in within_age_group:
            ij_data = infoes[((infoes["age1"] == i) & (infoes["age2"] == j)) | ((infoes["age1"] == j) & (infoes["age2"] == i))]
            edgeLen = len(ij_data)
            ij_data_gtex = new_gtex[((new_gtex["Age1"] == i) & (new_gtex["Age2"] == j)) | ((new_gtex["Age1"] == j) & (new_gtex["Age2"] == i))]
            all_edgeLen = len(ij_data_gtex)
            mean_int = round(edgeLen / all_edgeLen, 5) if all_edgeLen != 0 else 0
            labels.append(label)
            time_gap = within_age_group.get(f"{i}_VS_{j}", within_age_group.get(f"{j}_VS_{i}"))
            time_gaps.append(time_gap)
            mean_ints.append(mean_int)
        else:
            ij_data = infoes[((infoes["age1"] == i) & (infoes["age2"] == j)) | ((infoes["age1"] == j) & (infoes["age2"] == i))]
            edgeLen = len(ij_data)
            ij_data_gtex = new_gtex[((new_gtex["Age1"] == i) & (new_gtex["Age2"] == j)) | ((new_gtex["Age1"] == j) & (new_gtex["Age2"] == i))]
            all_edgeLen = len(ij_data_gtex)
            mean_int = round(edgeLen / all_edgeLen, 5) if all_edgeLen != 0 else 0
            labels.append(label)
            time_gap = abs(gene2time[i] - gene2time[j])
            time_gaps.append(time_gap)
            mean_ints.append(mean_int)

# Log-transform time gaps
log_time_gaps = np.log(time_gaps)

# Plotting
plt.figure()

colors = plt.cm.viridis(np.linspace(0, 1, 10))
colors = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"]

# Plot each point with unique label and color
for log_time_gap, mean_int, label, color in zip(log_time_gaps, mean_ints, labels, colors):
    plt.scatter(log_time_gap, mean_int, color=color, alpha=0.7, label=label)

plt.xlabel('Log(Time gaps)')
plt.ylabel('Loss percentage')

# Fit and plot regression line
m, b = np.polyfit(log_time_gaps, mean_ints, 1)
plt.plot(log_time_gaps, m * np.array(log_time_gaps) + b, color='red')

# Calculate Pearson correlation coefficient
corr, p_value = pearsonr(log_time_gaps, mean_ints)
print(corr, p_value)
plt.text(max(log_time_gaps) * 0.7, max(mean_ints) * 0.9, f'R = {corr:.2f}')

plt.legend()
plt.tight_layout()
plt.savefig("/data/chuand/essential_gene_evo/result/log_time_gap_vs_mean_interaction_loss2.pdf")

