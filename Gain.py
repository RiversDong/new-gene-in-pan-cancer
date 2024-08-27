import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import networkx as nx
from scipy.stats import pearsonr
import matplotlib as mpl

mpl.rcParams["pdf.fonttype"] = 42

def get_age(data, age_type):
    gene_age = data[["Gene", age_type]]
    gene_age = gene_age.dropna(axis=0, how="any")
    return gene_age

# File paths
tcga = "/data/chuand/essential_gene_evo/data/TCGA_interaction_type_top_tcga_mutual_20"
gtex = "/data/chuand/essential_gene_evo/data/GTEx_interaction_type_top_gtex_mutual_20"
age = "/data/chuand/essential_gene_evo/data/13059_2022_2821_MOESM2_ESM_Ma_gb_2022.txt"

# Load data
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
new_tcga = add_age(tcga, gene2age)

# Interaction data
out_res = "/data/chuand/essential_gene_evo/result/gained.csv"
age_groups = ["UC", "EM", "MM", "PSG"]
infoes = pd.read_csv(out_res, sep=",")
gene2time = {"PSG": 73, "MM": 180, "EM": 758, "UC": 4250}
time_gaps = []
mean_ints = []
labels = []

# For same-origin times, use mid-point value
within_age_group = {
    "UC_VS_UC": 1746,
    "EM_VS_EM": 289,
    "MM_VS_MM": 53.5,
    "PSG_VS_PSG": 36.5
}

# Dictionary to track unique pairs
pair_seen = {}

for i in age_groups:
    for j in age_groups:
        if f"{i}_VS_{j}" in pair_seen or f"{j}_VS_{i}" in pair_seen:
            continue  # Skip if the reverse pair has already been processed
        
        label = f"{i}-{j}({j}-{i})" if i < j else f"{j}-{i}"  # Unified label for symmetric pairs
        pair_seen[f"{i}_VS_{j}"] = True
        pair_seen[f"{j}_VS_{i}"] = True
        
        if f"{i}_VS_{j}" in within_age_group or f"{j}_VS_{i}" in within_age_group:
            ij_data = infoes[((infoes["age1"] == i) & (infoes["age2"] == j)) | ((infoes["age1"] == j) & (infoes["age2"] == i))]
            edgeLen = len(ij_data)
            ij_data_gtex = new_gtex[((new_gtex["Age1"] == i) & (new_gtex["Age2"] == j)) | ((new_gtex["Age1"] == j) & (new_gtex["Age2"] == i))]
            all_edgeLen_gtex = len(ij_data_gtex)
            mean_int = round(edgeLen / all_edgeLen_gtex, 5) if all_edgeLen_gtex != 0 else 0
            labels.append(label)
            time_gap = within_age_group.get(f"{i}_VS_{j}", within_age_group.get(f"{j}_VS_{i}"))
            time_gaps.append(time_gap)
            mean_ints.append(mean_int)
        else:
            ij_data = infoes[((infoes["age1"] == i) & (infoes["age2"] == j)) | ((infoes["age1"] == j) & (infoes["age2"] == i))]
            edgeLen = len(ij_data)
            ij_data_gtex = new_gtex[((new_gtex["Age1"] == i) & (new_gtex["Age2"] == j)) | ((new_gtex["Age1"] == j) & (new_gtex["Age2"] == i))]
            all_edgeLen_gtex = len(ij_data_gtex)
            mean_int = round(edgeLen / all_edgeLen_gtex, 5) if all_edgeLen_gtex != 0 else 0
            labels.append(label)
            time_gap = abs(gene2time[i] - gene2time[j])
            time_gaps.append(time_gap)
            mean_ints.append(mean_int)

# Take the logarithm of time_gaps
log_time_gaps = np.log(time_gaps)

# Plot log(time_gaps) vs mean_ints
plt.figure()

colors = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"]

# Scatter plot with transparency
for log_time_gap, mean_int, label, color in zip(log_time_gaps, mean_ints, labels, colors):
    plt.scatter(log_time_gap, mean_int, color=color, alpha=0.7, label=label)

plt.xlabel('Log (Time gaps)')
plt.ylabel('Gain percentage')

# Fit and plot a regression line
m, b = np.polyfit(log_time_gaps, mean_ints, 1)
plt.plot(log_time_gaps, m * np.array(log_time_gaps) + b, color='red')

# Calculate Pearson correlation coefficient
corr, p_value = pearsonr(log_time_gaps, mean_ints)
print(f"Pearson correlation: {corr}, P-value: {p_value}")
plt.text(max(log_time_gaps) * 0.7, max(mean_ints) * 0.9, f'R = {corr:.2f}')

plt.legend()
plt.tight_layout()
plt.savefig("/data/chuand/essential_gene_evo/result/time_gap_vs_mean_interaction_gain2.pdf")
