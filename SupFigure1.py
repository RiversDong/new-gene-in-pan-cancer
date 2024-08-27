import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ranksums
import seaborn as sns
import matplotlib
from numpy import mean
from scipy.stats import pearsonr

plt.rcParams.update({'font.family':'arial', 'font.size':13})

def get_age(data, age_type):
    gene_age = data[["Gene", age_type]]
    gene_age = gene_age
    gene_age = gene_age.dropna(axis=0, how="any")
    gene = list(gene_age["Gene"]); br = list(gene_age[age_type])
    gene2age = {}
    for i,j in zip(gene, br):
        gene2age[i] = j
    return gene2age

age = r"D:\essential_gene_evo\data\13059_2022_2821_MOESM2_ESM_Ma_gb_2022.txt"
findex = r'D:\essential_gene_evo\result\Sup.xlsx'

age_group = {"Cellular_organisms":1, "Eukaryota":2, "Opisthokonta":3, "Metazoa":4, "Eumetazoa":5, "Bilateria":6, "Chordata":7, "Euteleostomi":8, "Tetrapoda":9, "Amniota":10, "Mammalia":11, "Theria":12, "Eutheria":13, "Primate":14}
data = pd.read_csv(age, sep="\t", low_memory=False)
gene2age_14 = get_age(data, "14 age groups")
gene2age_4 = get_age(data, "4 merged age groups")

id_gene = data[["Gene", "Symbol"]]
ids = list(id_gene.iloc[1:]["Gene"].astype("str"))
genes = list(id_gene.iloc[1:]["Symbol"])
symbol2gene = {}
for i,j in zip(ids, genes):
    symbol2gene[j.upper()]=i

findex = pd.read_excel(findex, sheet_name="pan-score", names=["Gene", "pan-score"])
genes = findex["Gene"]; pan_scores = findex["pan-score"]
gene2pan_score = {}; ages = []; scores = []
age2boxplot={}
for i,j in zip(genes, pan_scores):
    pan_score = j
    i = i.upper()
    if i  in symbol2gene:
        iid = symbol2gene[i]
        if iid in gene2age_14:
            age = gene2age_14[iid]
            iid_age_group = age_group[age]
            ages.append(iid_age_group); scores.append(pan_score)
            if iid_age_group not in age2boxplot:
                age2boxplot[iid_age_group] = [pan_score]
            else:
                age2boxplot[iid_age_group].append(pan_score)
            #print(iid_age_group, pan_score)
data_list = []
branch = range(1,15)
data_list_four = []
group_1 = []; group_2=[]; group_3=[]; group_4=[]
for  i in branch:
    data_list.append(age2boxplot[i])
    if i in [1,2,3]:
        group_1.extend(age2boxplot[i])
    if i in [4,5,6,7,8,9,10]:
        group_2.extend(age2boxplot[i])
    if i in [11,12,13]:
        group_3.extend(age2boxplot[i])
    if  i in [14]:
        group_4.extend(age2boxplot[i])
data_list_four = [group_1, group_2, group_3, group_4]
fig, axes = plt.subplots(1, 1, figsize=(5.4, 4.8))

def ratio(group, cutoff):
    over = 0
    for i in group:
        if i >= cutoff:
            over = over + 1
    ratio = round(over/len(group), 3)
    return ratio
cutoff_list = [0.4, 0.5, 0.6, 0.7]

fig = plt.figure(figsize=(12, 4))
i_index = 1
x = list(range(1, 15))
feature=['Cellular_organisms', 'Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria', 'Chordata', 'Euteleostomi', 'Tetrapoda', 'Amniota', 'Mammalia', 'Theria', 'Eutheria', 'Primate']
for c in cutoff_list:
    ratios = []
    ax = fig.add_subplot(1,4,i_index)
    for i in data_list:
        i_ratio = ratio(i, c)
        ratios.append(i_ratio)
    markerline, stemlines, baseline = ax.stem(x, ratios, linefmt='g')
    plt.setp(stemlines, 'linewidth', 5)
    plt.setp(stemlines, 'color',["purple"])
    plt.setp(markerline, markersize = 6)
    plt.setp(baseline, visible=False)
    markerline.set_markerfacecolor('g')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_yticklabels("")
    max1 = max(ratios)
    max_index = ratios.index(max1)
    ax.set_yticks([max1+0.05])
    ax.grid(axis='y', linewidth=1.5, linestyle='dotted', color='purple')
    ax.set_xticks(x)
    ax.set_xticklabels(feature, rotation = 90)
    i_index = i_index +1
    print(c, ratios)
    plt.ylabel("Proporation")

plt.tight_layout()
#fig.savefig(r'D:\essential_gene_evo\new_figure_1\figure_s1_pan_cutoff_age_more_details.pdf', dpi=600, transparent=True)